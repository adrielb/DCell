#include "ImmersedInterfaceMethod.h"
#include "IIM_private.h"

void RemoveDuplicates( int *len, IIMIrregularNode *nodes[] ) ;

#undef __FUNCT__
#define __FUNCT__ "IIMUpdateSurfaceDerivatives_2D"
PetscErrorCode IIMUpdateSurfaceDerivatives_2D( IIM iim, LevelSet ls )
{
  int len, j;
  IIMIrregularNode *n, *nodes[iim->Np];
  PetscReal *eta, *xi;    // Local Coordinates in iim->lc
  PetscReal *s, *g;       // Arc length and surface quantity in lsq
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  ierr = PetscLogEventBegin(EVENT_IIMUpdateSurfaceDerivatives,0,0,0,0); CHKERRQ(ierr);
  LeastSqGetVecs(iim->lsq, &s, PETSC_NULL, &g, PETSC_NULL);
  
  int i;
  for( i = 0; i < ArrayLength(ls->irregularNodes); i++ )
  {
    ierr = ArrayGet(ls->irregularNodes,i,&n); CHKERRQ(ierr);
    ierr = SpatialIndexQueryPoints(iim->sidx, n->X, iim->eps, iim->Np, &len, (void**)nodes); CHKERRQ(ierr);
    RemoveDuplicates( &len, nodes);
    n->numNei = len;

    if( len < 3 ) {
      n->f1    = 0;
      n->f1_n  = 0;
      n->f1_nn = 0;
      n->f2    = 0;
      n->f2_n  = 0;
      n->f2_nn = 0;
      continue;
    }

    LocalCoorSetLength(iim->lc, len);
    LocalCoorSolve( iim->lc, iim->df, n, nodes);
    ierr = LeastSqSetNumPoints(iim->lsq, len); CHKERRQ(ierr);
    LocalCoorGetVecs(iim->lc, &eta, &xi, 0);
    for( j = 0; j < len; j++ )
    {
      s[j] = eta[j];
      g[j+0*len] = nodes[j]->F1;
      g[j+1*len] = nodes[j]->F2;
    }
    LeastSqSolve(iim->lsq);
    n->f1    = g[0+0*len];
    n->f1_n  = g[1+0*len];
    n->f1_nn = g[2+0*len];
    n->f2    = g[0+1*len];
    n->f2_n  = g[1+1*len];
    n->f2_nn = g[2+1*len];
  }
  
  ierr = PetscLogEventEnd(EVENT_IIMUpdateSurfaceDerivatives,0,0,0,0); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "IIMUpdateSurfaceDerivatives_3D"
PetscErrorCode IIMUpdateSurfaceDerivatives_3D( IIM iim, LevelSet ls )
{
  int len, i, j;
  IIMIrregularNode *n, *nodes[iim->Np];
  PetscReal *f1, *f2, *f3;
  PetscReal *s,  *r,  *g;
  PetscReal *ss, *rr, *nn;
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  ierr = PetscLogEventBegin(EVENT_IIMUpdateSurfaceDerivatives,0,0,0,0); CHKERRQ(ierr);
  LeastSqGetVecs(iim->lsq, &s, &r, &g, PETSC_NULL);
  
  for( i = 0; i < ArrayLength(ls->irregularNodes); i++ )
  {
    ierr = ArrayGet(ls->irregularNodes,i,(void*)&n); CHKERRQ(ierr);

    f1 = &n->f1;
    f2 = &n->f2;
    f3 = &n->f3;
    
    ierr = SpatialIndexQueryPoints(iim->sidx, n->X, iim->eps, iim->Np, &len, (void**)nodes); CHKERRQ(ierr);
    RemoveDuplicates( &len, nodes);
    n->numNei = len;
    if( len < 6 ) {
      for (j = 0; j < 6; ++j) {
        f1[j] = 0;
        f2[j] = 0;
        f3[j] = 0;
      }
      continue;
    }

    LocalCoorSetLength( iim->lc, len);
    LocalCoorSolve( iim->lc, iim->df, n, nodes);
    ierr = LeastSqSetNumPoints(iim->lsq, len); CHKERRQ(ierr);
    LocalCoorGetVecs( iim->lc, &ss, &nn, &rr );
    for( j = 0; j < len; j++ )
    {
      s[j] = ss[j];
      r[j] = rr[j];
      g[j+0*len] = nodes[j]->f1;
      g[j+1*len] = nodes[j]->f2;
      g[j+2*len] = nodes[j]->f3;
    }
    
    ierr = LeastSqSolve(iim->lsq); CHKERRQ(ierr);

    for( j = 0; j < 6; j++ ) {
      f1[j] = g[j+0*len];
      f2[j] = g[j+1*len];
      f3[j] = g[j+2*len];
    }
  }
  
  ierr = PetscLogEventEnd(EVENT_IIMUpdateSurfaceDerivatives,0,0,0,0); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

void RemoveDuplicates( int *len, IIMIrregularNode *nodes[] ) {
  int i,j;
  Coor d;
  const PetscReal tol = 1e-4;
  for (j = 0; j < *len; ++j) {
    for (i = j+1; i < *len; ++i) {
      d.x = nodes[i]->X.x - nodes[j]->X.x;
      d.y = nodes[i]->X.y - nodes[j]->X.y;
      d.z = nodes[i]->X.z - nodes[j]->X.z;
      if( d.x*d.x + d.y*d.y + d.z*d.z < tol ) {
        (*len)--;
        nodes[i] = nodes[*len];
        i--;
        continue;
      }
    }
  }
}
