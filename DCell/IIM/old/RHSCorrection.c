//#include "mkl.h"
#include "ImmersedInterfaceMethod.h"
#include "LevelSetMethod.h"

#undef __FUNCT__
#define __FUNCT__ "IIMComputeCorrection2D"
PetscInt EVENT_IIMComputeCorrection2D;
PetscErrorCode IIMComputeCorrection2D( IIM iim, JumpCondition Jump,
    LevelSet2D ls, Grid2D rhsC )
{
  int i,j;
  PetscReal a2,a4,a6,a8,a10,a12;
  PetscReal **phi = ls->g2d->v2;
  PetscInt len; //Number of nodes iim->kdtree returns 
  PetscReal *eta, *xi;    // Local Coordinates in iim->lc
  PetscReal *s, *g; // Arc param, Nodal Values and Surface Derivs in iim->lsq
  PetscReal w,dw,ddw,v,dv;
  IrregularNode *n, *node;
  GSList *slist, *iter;
  
  PetscErrorCode ierr;
  SETERRQ(PETSC_ERR_SUP, "IIMComputeCorrection2D deprecated, use IIM Simple");
  PetscFunctionBegin;
  PetscLogEventBegin(EVENT_IIMComputeCorrection2D,0,0,0,0);
//  PetscLogEventRegister(&EVENT_IIMComputeCorrection2D,"IIMComputeCorrection2D", 0);
  
  IIMUpdateIrregularNodes( iim , ls);
  LocalCoor2DGetVecs(iim->lc, &eta, &xi);
  LeastSqGetVecs(iim->lsq, &s, &g, PETSC_NULL);
  
  for( i = 0; i < ls->irregularNodes->len; i++ )
  {
    n = &g_array_index( ls->irregularNodes, IrregularNode, i);
    //Rotate the stencil into local coor space (eta, xi)
    LocalCoor2DSolveStencil( iim->lc, n );
    
    //calc a2 - a12 from 3.17
    a2=a4=a6=a8=a10=a12=0;
    for( j = 0; j < 5; j++)
    {
      if( phi[Stencil2Dx[j]+n->x][Stencil2Dy[j]+n->y] < 0)
        continue; // determine which points in the stencil are inside the interface
      a2 += DiscretePoisson[j];
      a4 += DiscretePoisson[j] * eta[j];
      a6 += DiscretePoisson[j] *  xi[j];
      a8 += DiscretePoisson[j] * eta[j] * eta[j];
      a10+= DiscretePoisson[j] *  xi[j] *  xi[j];
      a12+= DiscretePoisson[j] * eta[j] *  xi[j];
    }
    a8  /= two;
    a10 /= two;
    
    // reuse eta and xi, now global coordinates of ortho projections
    KDTreeRange(iim->kdtree, n, &slist );
    iter = slist;
    for( j = 0; j < iim->Np; j++ )
    {
      if( !iter ) break;
      node = KDTreeGetIrregularNode(iter);
      eta[j] = node->x + node->ox;
      xi[j]  = node->y + node->oy;
      iter = g_slist_next(iter);
    }
    len = j;
    ierr = LeastSqSetNumPoints(iim->lsq, len); CHKERRQ(ierr);

    // redefine eta and xi to local coordiates
    LocalCoor2DSetLength(iim->lc, len);
    LocalCoor2DSolve( iim->lc, n);
    LocalCoor2DToArcLength(iim->lc, s);
    
    Jump( iim->lsq, n, slist, &w, &dw, &ddw, &v, &dv );
    
    n->c = a2 * w + a12 * dv + ( a6 + a12 * n->k ) * dw + a10 * ddw +
            v * ( a4 + (a8-a10) * n->k ) - a8 * ddw;
    
    rhsC->v2[n->x][n->y] += n->c;
  }
  PetscLogEventEnd(EVENT_IIMComputeCorrection2D,0,0,0,0);
  PetscFunctionReturn(0);
}

void IIMDiscreteCompatabilityCondition( LevelSet2D ls, Grid2D c)
{
  IrregularNode *n;
  PetscReal avg;
  VecSum(c->v,&avg);
  avg /= ls->irregularNodes->len;
  for( int i = 0; i < ls->irregularNodes->len; i++ )
  {
    n = &g_array_index( ls->irregularNodes, IrregularNode, i);
    c->v2[n->x][n->y] -= avg;
  }
}