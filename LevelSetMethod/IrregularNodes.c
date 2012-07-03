#include "LevelSetMethod.h"
#include "LSM_private.h"

/*
 * NOTE: An irregular grid point cannot be on the boundary of the grid!!!
 */
#undef __FUNCT__
#define __FUNCT__ "LevelSetUpdateIrregularNodeList"
PetscErrorCode LevelSetUpdateIrregularNodeList( LevelSet ls )
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscLogEventBegin(EVENT_LevelSetUpdateIrregularNodeList,0,0,0,0); CHKERRQ(ierr);
  
  ierr = ArraySetSize(ls->irregularNodes, 0); CHKERRQ(ierr);
  ierr = ArrayZero(ls->irregularNodes); CHKERRQ(ierr);

  if( ls->phi->is2D )
  {
    ierr = LevelSetUpdateIrregularNodeList_2D( ls ); CHKERRQ(ierr);
  } else {
    ierr = LevelSetUpdateIrregularNodeList_3D( ls ); CHKERRQ(ierr);
  }
  
  ierr = PetscLogEventEnd(EVENT_LevelSetUpdateIrregularNodeList,0,0,0,0); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "LevelSetUpdateIrregularNodeList_2D"
PetscErrorCode LevelSetUpdateIrregularNodeList_2D( LevelSet ls )
{
  int i, j, k, I, J, ni, nj;
  const int numNei = 4;
  const int nei[][2] = {{1,0},{-1,0},{0,-1},{0,1}};
  IrregularNode *n;
  PetscReal **phi, phiHI, phiLO;
  PetscReal sten[3][3];
  PetscReal local[5][5];
  iCoor *band;
  int b;
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  ierr = GridGet(ls->phi,&phi); CHKERRQ(ierr);
  for( b = 0; b < ArrayLength(ls->band); ++b)
  {
    ierr = ArrayGet(ls->band,b,&band); CHKERRQ(ierr);
    i = band->x;
    j = band->y;
    
    // Cell-centered Irregular Node
    for( J = -2; J <= 2; ++J) {
      for( I = -2; I <= 2; ++I) {
        local[J+2][I+2] = phi[J+j][I+i];
      }
    } // for local 5x5 stencil

    // Cell-centered Irregular Node
    for( J = -1; J < 2; ++J)
    {
      for( I = -1; I < 2; ++I)
      {
        sten[J+1][I+1] = phi[J+j][I+i];
      }
    } // for local 3x3 stencil

    // Add ortho-proj for FMM boundary condition
    for( k = 0; k < numNei; ++k)
    {
      ni = 1 + nei[k][0];
      nj = 1 + nei[k][1];
      if( sten[1][1] * sten[nj][ni] <= 0. )
      {
        ierr = ArrayAppend( ls->irregularNodes, &n ); CHKERRQ(ierr);
        OrthogonalProjection2D( sten, local, &n->op);
        n->pos = *band;
        n->axis  = -1; // no-axis
        n->shift = -1;
        n->signCenter = sten[1][1] > 0. ? 1 : -1;
        n->X.x = n->pos.x + n->op.x;
        n->X.y = n->pos.y + n->op.y;
        break;
      }
    }

    // Add IIM irregular grid point
    // Cell-centered gradients
    for( k = U_FACE; k <= V_FACE; ++k)
    {
      ni = 1 + STAGGERED_GRID[k].x;
      nj = 1 + STAGGERED_GRID[k].y;
      if( sten[1][1] * sten[nj][ni] <= 0. )
      {
        ierr = ArrayAppend( ls->irregularNodes, &n ); CHKERRQ(ierr);
        n->d = sten[1][1] / (sten[1][1] - sten[nj][ni]);
        n->signCenter = sten[1][1] > 0. ? 1 : -1;
        n->signFace = n->d < 0.5 ? -n->signCenter : n->signCenter;
        n->pos = *band;
        n->shift = CELL_CENTER;
        n->axis  = k-U_FACE; // assuming U_FACE == 1, x-axis is 0, y-axis is 1
        n->X.x = n->pos.x + n->d * STAGGERED_GRID[k].x;
        n->X.y = n->pos.y + n->d * STAGGERED_GRID[k].y;
      } // if irreg
    } // for k in {U_FACE,V_FACE}

    // U-Velocity Laplacian
    phiHI = ( sten[1][0] + sten[1][1] ) / 2.;
    phiLO = ( sten[0][0] + sten[0][1] ) / 2.;
    if( phiHI * phiLO <= 0. )
    {
      ierr = ArrayAppend( ls->irregularNodes, &n ); CHKERRQ(ierr);
      n->d = phiHI / (phiHI - phiLO);
      n->signCenter = phiHI > 0. ? 1 : -1;
      n->signFace = n->d < 0.5 ? -n->signCenter : n->signCenter;
      n->pos = *band;
      n->shift = U_FACE;
      n->axis  = V_FACE-U_FACE; // y-axis == 1
      n->X.x = n->pos.x + STAGGERED_GRID[U_FACE].x / 2.;
      n->X.y = n->pos.y + STAGGERED_GRID[U_FACE].y / 2. - n->d;
    } // if irreg

    // V-Velicty Laplacian
    phiHI = ( sten[0][1] + sten[1][1] ) / 2.;
    phiLO = ( sten[0][0] + sten[1][0] ) / 2.;
    if( phiHI * phiLO <= 0. )
    {
      ierr = ArrayAppend( ls->irregularNodes, &n ); CHKERRQ(ierr);
      n->d = phiHI / (phiHI - phiLO);
      n->signCenter = phiHI > 0. ? 1 : -1;
      n->signFace = n->d < 0.5 ? -n->signCenter : n->signCenter;
      n->pos = *band;
      n->shift = V_FACE;
      n->axis  = U_FACE-U_FACE; // x-axis == 0
      n->X.x = n->pos.x + STAGGERED_GRID[V_FACE].x / 2. - n->d;
      n->X.y = n->pos.y + STAGGERED_GRID[V_FACE].y / 2.;
    } // if irreg
  } // for b in band
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "LevelSetWriteIrregularNodeList_2D"
PetscErrorCode LevelSetWriteIrregularNodeList_2D( Array irregularNodes, PetscViewer viewer )
{
  int i,j;
  int len = ArrayLength( irregularNodes );
  const int rowlen = 10;
  PetscReal row[rowlen]; // { X, Y, nv, nv, f1, f2, k, f1_n }
  IrregularNode *node, *nodes = ArrayGetData( irregularNodes );
  PetscErrorCode ierr;

  PetscFunctionBegin;
  for ( i = 0; i < len; ++i) {
    node = &nodes[i];
    if( node->axis == -1 ) continue;
//    if( (node->shift == -1 && node->axis == -1) )
    {
      j=0;
      row[j++] = node->X.x;
      row[j++] = node->X.y;
      row[j++] = node->pos.x;
      row[j++] = node->pos.y;
      row[j++] = node->nx;
      row[j++] = node->ny;
      row[j++] = node->f1;
      row[j++] = node->fa1;
      row[j++] = node->k;
      row[j++] = node->f1_n;
      ierr = PetscViewerBinaryWrite(viewer,row,rowlen,PETSC_DOUBLE,PETSC_TRUE); CHKERRQ(ierr);
    }
  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "LevelSetWriteIrregularNodeList_3D"
PetscErrorCode LevelSetWriteIrregularNodeList_3D( Array irregularNodes, PetscViewer viewer )
{
  int i,j;
  int len = ArrayLength( irregularNodes );
  const int rowlen = 23;
  PetscReal row[rowlen]; // { X, Y, Z, nv, nv, nv, k, f1 }
  IrregularNode *node, *nodes = ArrayGetData( irregularNodes );
  PetscErrorCode ierr;

  PetscFunctionBegin;
  for ( i = 0; i < len; ++i) {
    node = &nodes[i];
//    if( (node->shift == -1 && node->axis == -1) )
    {
      j=0;
      row[j++] = node->X.x;    // 0
      row[j++] = node->X.y;
      row[j++] = node->X.z;
      row[j++] = node->pos.x;  // 3
      row[j++] = node->pos.y;
      row[j++] = node->pos.z;
      row[j++] = node->nx;     // 6
      row[j++] = node->ny;
      row[j++] = node->nz;
      row[j++] = node->k;      // 9
      row[j++] = node->f1;
      row[j++] = node->f1_n;
      row[j++] = node->f1_t;   //12
      row[j++] = node->f1_nn;
      row[j++] = node->f1_tt;
      row[j++] = node->f1_nt;  //15
      row[j++] = node->sx;     //16
      row[j++] = node->sy;
      row[j++] = node->sz;
      row[j++] = node->rx;     //19
      row[j++] = node->ry;
      row[j++] = node->rz;     //21
      row[j++] = node->numNei; //22
      ierr = PetscViewerBinaryWrite(viewer,row,rowlen,PETSC_DOUBLE,PETSC_TRUE); CHKERRQ(ierr);
    }
  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "LevelSetWriteIrregularNodeList"
PetscErrorCode LevelSetWriteIrregularNodeList( LevelSet ls, int idx )
{
  char filename[PETSC_MAX_PATH_LEN], wd[PETSC_MAX_PATH_LEN];
  PetscViewer viewer;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscLogEventBegin(EVENT_LevelSetWriteIrregularNodeList,0,0,0,0); CHKERRQ(ierr);

  // TODO: label with levelset id tag:
//  ierr = ArrayWrite( irreg, "irregNode", idx ); CHKERRQ(ierr);

  ierr = DCellGetWorkingDirectory(wd); CHKERRQ(ierr);
  ierr = PetscSNPrintf(filename,PETSC_MAX_PATH_LEN,"%s/%s.irregNode.%d.array",wd,ls->phi->name,idx+FILE_COUNT_START); CHKERRQ(ierr);
  ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF,filename,FILE_MODE_WRITE,&viewer); CHKERRQ(ierr);
  if( ls->phi->is2D ) {
    ierr = LevelSetWriteIrregularNodeList_2D( ls->irregularNodes, viewer ); CHKERRQ(ierr);
  } else {
    ierr = LevelSetWriteIrregularNodeList_3D( ls->irregularNodes, viewer ); CHKERRQ(ierr);
  }
  ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);
  ierr = PetscLogEventEnd(EVENT_LevelSetWriteIrregularNodeList,0,0,0,0); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
