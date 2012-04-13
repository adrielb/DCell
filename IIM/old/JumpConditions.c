#include "ImmersedInterfaceMethod.h"

/* Pressure Jump Condition
 * [p] = f1
 * [pn]= df2/ds
 */
void JumpConditionPressure( LeastSq lsq, IrregularNode *n, GSList *slist,
    PetscReal *w, PetscReal *dw, PetscReal *ddw, PetscReal *v, PetscReal *dv)
{
  PetscReal *g, len;
  GSList *iter;
  
  LeastSqGetVecs(lsq, PETSC_NULL, &g, &len);
  
  // Calculate and cache normal force
  iter = slist;
  for( int i = 0; i < len; i++ )
  {
    g[i] = KDTreeGetIrregularNode(iter)->f1;
    iter = g_slist_next(iter);
  }
  LeastSqSolve(lsq);
//  n->f1   = g[0];
  n->df1  = g[1];
  n->ddf1 = g[2];
  
  // Calculate and cache tangential force
  iter = slist;
  for( int i = 0; i < len; i++ )
  {
    g[i] = KDTreeGetIrregularNode(iter)->f2;
    iter = g_slist_next(iter);
  }
  LeastSqSolve(lsq);
//  n->f2   = g[0];
  n->df2  = g[1];
  n->ddf2 = g[2];
  
  *w   = n->f1;
  *dw  = n->df1;
  *ddw = n->ddf1;
  *v   = n->df2;
  *dv  = n->ddf2;
}

/* Velocity Jump Condition
 * [u] = 0
 * [un]= f2 sin(q)
 */
void JumpConditionXVelocity( LeastSq lsq, IrregularNode *n, GSList *slist,
    PetscReal *w, PetscReal *dw, PetscReal *ddw, PetscReal *v, PetscReal *dv)
{
  PetscReal *g, len;
  IrregularNode *node;
  GSList *iter=slist;
  
  *w=*dw=*ddw=0;
  
  LeastSqGetVecs(lsq, PETSC_NULL, &g, &len);
  
  for( int i = 0; i < len; i++ )
  {
    node = KDTreeGetIrregularNode(iter);
    g[i] = node->f2 * node->ny;
    iter = g_slist_next(iter);
  }
  
  LeastSqSolve(lsq);

  *v= n->f2 * n->ny;
  *dv= g[1];
}

/* Velocity Jump Condition
 * [v] = 0
 * [vn]= f2 cos(q)
 */
void JumpConditionYVelocity( LeastSq lsq, IrregularNode *n, GSList *slist,
    PetscReal *w, PetscReal *dw, PetscReal *ddw, PetscReal *v, PetscReal *dv)
{
  PetscReal *g, len;
  IrregularNode *node;
  GSList *iter=slist;
  
  *w=*dw=*ddw=0;
  
  LeastSqGetVecs(lsq, PETSC_NULL, &g, &len);
  
  for( int i = 0; i < len; i++ )
  {
    node = KDTreeGetIrregularNode(iter);
    g[i] = -1 * node->f2 * node->nx;
    iter = g_slist_next(iter);
  }
  
  LeastSqSolve(lsq);
  
  *v= -1 * n->f2 * n->nx;
  *dv= g[1];
}