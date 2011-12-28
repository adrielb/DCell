#include "ImmersedInterfaceMethod.h"
#include "IIM_private.h"

#undef __FUNCT__
#define __FUNCT__ "IIMRHSUpdate"
PetscErrorCode IIMUpdateRHS( IIM iim, LevelSet ls, int ga )
{
  int i;
  IrregularNode *node;
  int len = ArrayLength(ls->irregularNodes);
  Jump jump;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscLogEventBegin(EVENT_IIMUpdateRHS,0,0,0,0); CHKERRQ(ierr);
  ierr = IIMUpdateSurfaceQuantities(iim,ls); CHKERRQ(ierr);

  ierr = ArraySetSize(iim->coor,0); CHKERRQ(ierr);
  ierr = ArraySetSize(iim->idx,0); CHKERRQ(ierr);
  ierr = ArraySetSize(iim->val,0); CHKERRQ(ierr);

  for( i = 0; i < len; i++ )
  {
    ierr = ArrayGet(ls->irregularNodes,i,(void*)&node); CHKERRQ(ierr);
    JumpPressure( node, &jump );

    // C{ px = p[i] - p[i-1] }
    ierr = IIMPressureGradientCorrection( iim, node, jump ); CHKERRQ(ierr);

    if( node->shift == CELL_CENTER ) {
      JumpVelocity( *iim->mu, node, &jump, node->axis );
    } else {
      JumpVelocity( *iim->mu, node, &jump, node->shift - U_FACE );
    }

    // C{ laplace(u) }
    ierr = IIMLaplaceCorrection( iim, node, jump ); CHKERRQ(ierr);

    // C{ ux = u[i+1] - u[i] }
    ierr = IIMVelocityGradientCorrection( iim, node, jump ); CHKERRQ(ierr);
  }


  len = ArrayLength(iim->coor);
  ierr = ArraySetSize(iim->idx,len); CHKERRQ(ierr);
  int** idx = ArrayGetData(iim->idx);
  int*  coor;
  for (i = 0; i < len; ++i) {
    ierr = ArrayGet(iim->coor,i,&coor); CHKERRQ(ierr);
    idx[i] = coor;
//    printf("{%d:%d,%d,%d},",coor[0],coor[1],coor[2],coor[3]);
  }

  ierr = GAScatterAcc(ga,iim->idx,iim->val); CHKERRQ(ierr);
  ierr = PetscLogEventEnd(EVENT_IIMUpdateRHS,0,0,0,0); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
