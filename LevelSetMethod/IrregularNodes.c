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
  PetscReal **phi;
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
    for( J = -1; J <= 1; ++J)
    {
      for( I = -1; I <= 1; ++I)
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
        n->sign = sten[1][1] > 0. ? 1 : -1;
        n->X.x = n->pos.x + n->op.x;
        n->X.y = n->pos.y + n->op.y;
        break;
      }
    }

  } // for b in band
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "LevelSetWriteIrregularNodeList_2D"
PetscErrorCode LevelSetWriteIrregularNodeList_2D( Array irregularNodes, PetscViewer viewer )
{
  int i;
  int rowlen;
  const int len = ArrayLength( irregularNodes );
  IrregularNode *nodes = ArrayGetData( irregularNodes );
  IrregularNode *node;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  for ( i = 0; i < len; ++i) {
    node = &nodes[i];

    PetscReal row[] = {
        node->X.x,
        node->X.y,
        node->pos.x,
        node->pos.y
    };
    rowlen = sizeof(row)/sizeof(PetscReal);
    ierr = PetscViewerBinaryWrite(viewer,row,rowlen,PETSC_DOUBLE,PETSC_TRUE); CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "LevelSetWriteIrregularNodeList_3D"
PetscErrorCode LevelSetWriteIrregularNodeList_3D( Array irregularNodes, PetscViewer viewer )
{
  int i;
  int rowlen;
  const int len = ArrayLength( irregularNodes );
  IrregularNode *nodes = ArrayGetData( irregularNodes );
  IrregularNode *node;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  for ( i = 0; i < len; ++i) {
    node = &nodes[i];
    PetscReal row[] = {
        node->X.x,    // 0
        node->X.y,
        node->X.z,
        node->pos.x,  // 3
        node->pos.y,
        node->pos.z
    };
    rowlen = sizeof(row)/sizeof(PetscReal);
    ierr = PetscViewerBinaryWrite(viewer,row,rowlen,PETSC_DOUBLE,PETSC_TRUE); CHKERRQ(ierr);
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
