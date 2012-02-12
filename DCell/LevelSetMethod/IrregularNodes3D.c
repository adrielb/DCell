#include "LevelSetMethod.h"
#include "LSM_private.h"

#undef __FUNCT__
#define __FUNCT__ "LevelSetUpdateIrregularNodeList_3D"
PetscErrorCode LevelSetUpdateIrregularNodeList_3D( LevelSet ls )
{
  int i, j, k;
  int I, J, K;
  int ni, nj, nk;
  int m;
  const int numNei = 6;
  IrregularNode *n;
  PetscReal ***phi, phiHI, phiLO;
  PetscReal sten[3][3][3];
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
    k = band->z;

    // Cell-centered Irregular Node
    for( K = -1; K < 2; ++K)
    {
      for( J = -1; J < 2; ++J)
      {
        for( I = -1; I < 2; ++I)
        {
          sten[K+1][J+1][I+1] = phi[K+k][J+j][I+i];
        }
      }
    }// for local 3x3x3 stencil

    // Add ortho-proj for FMM boundary condition
    for( m = 0; m < numNei; ++m)
    {
      ni = 1 + STAR[m][0];
      nj = 1 + STAR[m][1];
      nk = 1 + STAR[m][2];
      if( sten[1][1][1] * sten[nk][nj][ni] <= 0. )
      {
        ierr = ArrayAppend( ls->irregularNodes, &n ); CHKERRQ(ierr);
        OrthogonalProjection3D( sten, &n->op);
        n->pos = *band;
        n->axis  = -1; // no-axis
        n->shift = -1;
        n->signCenter = sten[1][1][1] > 0. ? 1 : -1;
        n->X.x = n->pos.x + n->op.x;
        n->X.y = n->pos.y + n->op.y;
        n->X.z = n->pos.z + n->op.z;
        break;
      }
    }

    // Add IIM irregular grid point
    // Cell-centered gradients
    for( m = U_FACE; m <= W_FACE; ++m)
    {
      ni = 1 + STAGGERED_GRID[m].x;
      nj = 1 + STAGGERED_GRID[m].y;
      nk = 1 + STAGGERED_GRID[m].z;
      if( sten[1][1][1] * sten[nk][nj][ni] <= 0. )
      {
        ierr = ArrayAppend( ls->irregularNodes, &n ); CHKERRQ(ierr);
        n->d = sten[1][1][1] / (sten[1][1][1] - sten[nk][nj][ni]);
        n->signCenter = sten[1][1][1] > 0. ? 1 : -1;
        n->signFace = n->d < 0.5 ? -n->signCenter : n->signCenter;
        n->pos = *band;
        n->shift = CELL_CENTER;
        n->axis  = m-U_FACE; // assuming U_FACE == 1, x-axis is 0, y-axis is 1
        n->X.x = n->pos.x + n->d * STAGGERED_GRID[m].x;
        n->X.y = n->pos.y + n->d * STAGGERED_GRID[m].y;
        n->X.z = n->pos.z + n->d * STAGGERED_GRID[m].z;
      } // if irreg
    } // for m in {U_FACE,V_FACE,W_FACE}

    VelFace face, axis;
    for (face = U_FACE; face <= W_FACE; ++face) {
      iCoor g = STAGGERED_GRID[face];
      phiHI = ( sten[1+g.z][1+g.y][1+g.x] + sten[1][1][1] ) / 2.;
      for (axis = U_FACE; axis < W_FACE; ++axis) {
        if( face == axis ) continue;
        iCoor lo = STAGGERED_GRID[axis];
        phiLO = ( sten[1+g.z+lo.z][1+g.y+lo.y][1+g.x+lo.x] + sten[1+lo.z][1+lo.y][1+lo.x] ) / 2.;
        if( phiHI * phiLO <= 0. )
        {
          ierr = ArrayAppend( ls->irregularNodes, &n ); CHKERRQ(ierr);
          n->d = phiHI / (phiHI - phiLO);
          n->signCenter = phiHI > 0. ? 1 : -1;
          n->signFace = n->d < 0.5 ? -n->signCenter : n->signCenter;
          n->pos = *band;
          n->shift = face;
          n->axis  = axis - U_FACE;
          n->X.x = n->pos.x + STAGGERED_GRID[face].x / 2. - n->d * STAGGERED_GRID[axis].x;
          n->X.y = n->pos.y + STAGGERED_GRID[face].y / 2. - n->d * STAGGERED_GRID[axis].y;
          n->X.z = n->pos.z + STAGGERED_GRID[face].z / 2. - n->d * STAGGERED_GRID[axis].z;
        } // if irreg
      } // axis
    } // face
  } // for b in band
  PetscFunctionReturn(0);
}
