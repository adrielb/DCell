#include "ImmersedInterfaceMethod.h"
#include "IIM_private.h"

#undef __FUNCT__
#define __FUNCT__ "IIMUpdateSurfaceDerivatives_2D"
PetscErrorCode IIMUpdateSurfaceDerivatives_2D( IIM iim, LevelSet ls )
{
  IrregularNode *n, *nodes[iim->Np];
  int len, j;
  PetscReal *eta, *xi;    // Local Coordinates in iim->lc
  PetscReal *s, *g;       // Arc length and surface quantity in lsq
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  ierr = PetscLogEventBegin(EVENT_IIMUpdateSurfaceDerivatives,0,0,0,0); CHKERRQ(ierr);
  LocalCoor2DGetVecs(iim->lc, &eta, &xi);
  LeastSqGetVecs(iim->lsq, &s, PETSC_NULL, &g, PETSC_NULL);
  
  int i;
  for( i = 0; i < ArrayLength(ls->irregularNodes); i++ )
  {
    ierr = ArrayGet(ls->irregularNodes,i,&n); CHKERRQ(ierr);
    ierr = IIMSurfaceNeighbors_2D(iim,ls,n,(IrregularNode**)&nodes,&len); CHKERRQ(ierr);

    // TODO: this is a hack, need to ensure that all nodes are not islands ( have neighbors in IrregNodes update )
    if( len < 3 ) {
      n->f1    = 0;
      n->f1_n  = 0;
      n->f1_nn = 0;
      n->f2    = 0;
      n->f2_n  = 0;
      n->f2_nn = 0;
      continue;
    }

    for( j = 0; j < len; j++ )
    {
      eta[j] = nodes[j]->X.x;
      xi[j]  = nodes[j]->X.y;
    }
    
    LocalCoorSetLength(iim->lc, len);
    LocalCoor2DSolve( iim->lc, iim->dh, n);
    
    for( j = 0; j < len; j++ )
    {
//      LocalCoor2DToArcLength(iim->lc, n, nodes[j]->nx, nodes[j]->ny, j, &s[j]);
      s[j] = -xi[j];
    }
    
    ierr = LeastSqSetNumPoints(iim->lsq, len); CHKERRQ(ierr);

    for( j = 0; j < len; j++ )
    {
      g[j] = nodes[j]->f1;
    }
    LeastSqSolve(iim->lsq);
//    n->f1    = g[0];
    n->f1_n  = g[1];
    n->f1_nn = g[2];
    
    for( j = 0; j < len; j++ )
    {
      g[j] = nodes[j]->f2;
    }
    LeastSqSolve(iim->lsq);
//    n->f2    = g[0];
    n->f2_n  = g[1];
    n->f2_nn = g[2];
  }
  
  ierr = PetscLogEventEnd(EVENT_IIMUpdateSurfaceDerivatives,0,0,0,0); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
/*
 *
 *
 *
 *
#undef __FUNCT__
#define __FUNCT__ "IIMUpdateSurfaceDerivatives_3D"
PetscErrorCode IIMUpdateSurfaceDerivatives_3D( IIM iim, LevelSet ls )
{
  int len, j;
  IrregularNode *n, *nodes[iim->Np];
  PetscReal *f1, *f2, *f3;
  PetscReal *s, *r, *g;
  PetscReal *ee, *nn, *tt;
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  
  LeastSqGetVecs(iim->lsq, &s, &r, &g, PETSC_NULL);
  LocalCoor3DGetVecs( iim->lc, &nn, &tt, &ee );
  
  int i;
  for( i = 0; i < ArrayLength(ls->irregularNodes); i++ )
  {
    ierr = ArrayGet(ls->irregularNodes,i,(void*)&n); CHKERRQ(ierr);

    f1 = &n->f1;
    f2 = &n->f2;
    f3 = &n->f3;
    
    ierr = IIMSurfaceNeighbors_3D(iim,ls->g->n,ls->irregularNodeGrid,n,(IrregularNode**)&nodes,&len); CHKERRQ(ierr);

    for( j = 0; j < len; j++ )
    {
      ee[j] = nodes[j]->x + nodes[j]->ox;
      nn[j] = nodes[j]->y + nodes[j]->oy;
      tt[j] = nodes[j]->z + nodes[j]->oz;
    }
    
    LocalCoorSetLength(iim->lc, len);
    LocalCoor3DSolve( iim->lc, n);
    
    for( j = 0; j < len; j++ )
    {
      s[j] = nn[j];
      r[j] = tt[j];
    }
    
    ierr = LeastSqSetNumPoints(iim->lsq, len); CHKERRQ(ierr);

    for( j = 0; j < len; j++ )
      g[j] = nodes[j]->f1;
    LeastSqSolve(iim->lsq);
    for( j = 1; j < 5; j++ )
      f1[j] = g[j];
    
    for( j = 0; j < len; j++ )
      g[j] = nodes[j]->f2;
    LeastSqSolve(iim->lsq);
    for( j = 1; j < 5; j++ )
      f2[j] = g[j];
    
    for( j = 0; j < len; j++ )
      g[j] = nodes[j]->f3;
    LeastSqSolve(iim->lsq);
    for( j = 1; j < 5; j++ )
      f3[j] = g[j];
  }
  
  PetscFunctionReturn(0);
}
*/
