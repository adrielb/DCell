#include "DWorld.h"


#undef __FUNCT__
#define __FUNCT__ "RHSJacobian"
PetscInt EVENT_RHSJacobian;
PetscErrorCode RHSJacobian(TS ts,PetscReal t,Vec global_in,Mat *JJ,Mat *PP,
    MatStructure *flag,void *ptr)
{
  DWorld ctx = (DWorld)ptr;
  DA da = ctx->da;
  DALocalInfo i;
  Reaction rxn = ctx->rxnWorld;
  Mat J = *JJ;
  MatStencil row, col;
  PetscReal ****chem;
  PetscInt c;
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  PetscLogEventBegin(EVENT_RHSJacobian,0,0,0,0);

//  ierr = MatCopy(ctx->L,ctx->J,SAME_NONZERO_PATTERN); CHKERRQ(ierr);
  ierr = MatCopy(ctx->L,ctx->J,DIFFERENT_NONZERO_PATTERN); CHKERRQ(ierr);
  ierr = DAGetLocalInfo(da, &i); CHKERRQ(ierr);

  ierr = DAVecGetArrayDOF(ctx->da,ctx->global,&chem);CHKERRQ(ierr);  
  
  /* Update block diagonal part of Jacobian (the reaction part) */

  for (row.k = i.zs; row.k < i.zs+i.zm; ++row.k)
  {
    for (row.j = i.ys; row.j < i.ys+i.ym; ++row.j)
    {
      for (row.i = i.xs; row.i < i.xs+i.xm; ++row.i)
      {
       // if( rxntype[row.k][row.j][row.i] )
       
        ReactionUpdateJacobian(rxn, chem[row.k][row.j][row.i]);
        col.i = row.i;
        col.j = row.j;
        col.k = row.k;
        for( c = 0; c < rxn->jac_length; ++c)
        {
          row.c = rxn->rows[c];
          col.c = rxn->cols[c];
          MatSetValuesStencil(J,1,&row,1,&col,&rxn->jac[c],ADD_VALUES);
        }
//        MatSetValuesBlockedStencil(J, 1, &row, 1, &row, rxn->jac, ADD_VALUES );
      }
    }
  }
  /*
  for (s = 0; s < lsLen; ++s) // for each cell
  {
    for ( (i,j,k) local ls coor )
    {
      if( ls[s]->v3[k][j][i] < 0 // IsInside)
      {
        dcell[s]
      }
    }
    
    for (r = 0; r < ls[s]->irregularNodes->len; ++r) // for membrane rxns
    {
      n = get(r);
      
    }
  }
  */
  ierr = DAVecRestoreArrayDOF(da,ctx->global,&chem);CHKERRQ(ierr);
  
  ierr = MatAssemblyBegin(J,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(J,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  
//  *flag = DIFFERENT_NONZERO_PATTERN;
  *flag = SAME_NONZERO_PATTERN;
  
  PetscLogEventEnd(EVENT_RHSJacobian,0,0,0,0);
  PetscFunctionReturn(0);
}