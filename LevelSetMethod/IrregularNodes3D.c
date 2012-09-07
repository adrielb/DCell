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
  PetscReal ***phi;
  PetscReal  sten[3][3][3];
  PetscReal local[5][5][5];
  int b;
  const int len = ArrayLength(ls->band);
  const iCoor *band = ArrayGetData(ls->band);
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = GridGet(ls->phi,&phi); CHKERRQ(ierr);
  for( b = 0; b < len; ++b)
  {
    i = band[b].x;
    j = band[b].y;
    k = band[b].z;

    // Cell-centered Irregular Node
    for( K = -2; K <= 2; ++K) {
      for( J = -2; J <= 2; ++J) {
        for( I = -2; I <= 2; ++I) {
          local[K+2][J+2][I+2] = phi[K+k][J+j][I+i];
        }
      }
    } // for local 5x5x5 stencil

    // Cell-centered Irregular Node
    for( K = -1; K <= 1; ++K) {
      for( J = -1; J <= 1; ++J) {
        for( I = -1; I <= 1; ++I) {
          sten[K+1][J+1][I+1] = phi[K+k][J+j][I+i];
        }
      }
    } // for local 3x3x3 stencil

    // Add ortho-proj for FMM boundary condition
    for( m = 0; m < numNei; ++m)
    {
      ni = 1 + STAR[m][0];
      nj = 1 + STAR[m][1];
      nk = 1 + STAR[m][2];
      if( sten[1][1][1] * sten[nk][nj][ni] <= 0. )
      {
        ierr = ArrayAppend( ls->irregularNodes, &n ); CHKERRQ(ierr);
        OrthogonalProjection3D( sten, local, &n->op);
        n->pos = band[b];
        n->sign = sten[1][1][1] > 0. ? 1 : -1;
        n->X.x = n->pos.x + n->op.x;
        n->X.y = n->pos.y + n->op.y;
        n->X.z = n->pos.z + n->op.z;
        break;
      }
    }

  } // for b in band
  PetscFunctionReturn(0);
}
