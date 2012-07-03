#include "ImmersedInterfaceMethod.h"

void MyForce( IrregularNode *n, void *ctx );

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc, char **args)
{
  PetscErrorCode  ierr;
  PetscFunctionBegin;
  ierr = DCellInit(); CHKERRQ(ierr);

  PetscReal dx = 1;
  Coor dh = {dx, dx, 0};
  Jump jump = { .j = 0, .x = 0, .y = 0, .z = 0, .xx = 0, .yy = 0, .zz = 0 };

  LevelSet ls;
  ierr = LevelSetCreate(dh,(iCoor){-2,-2,0},(iCoor){6,6,0},&ls); CHKERRQ(ierr);
  PetscReal **phi;
  ierr = VecSet( ls->phi->v, 1 ); CHKERRQ(ierr);
  ierr = GridGet( ls->phi, &phi); CHKERRQ(ierr);
  phi[0][0] = -1;
  iCoor *band;
  ierr = ArrayAppend(ls->band,&band); CHKERRQ(ierr);
  band->x = 0;
  band->y = 0;
  ierr = ArrayAppend(ls->band,&band); CHKERRQ(ierr);
  band->x = 0;
  band->y = 1;
  ierr = ArrayAppend(ls->band,&band); CHKERRQ(ierr);
  band->x = 1;
  band->y = 0;
  ierr = ArrayAppend(ls->band,&band); CHKERRQ(ierr);
  band->x = 1;
  band->y = 1;

  IIM iim;
  ierr = IIMCreate(PETSC_TRUE,20,dh,&iim); CHKERRQ(ierr);
  ierr = IIMSetViscosity(iim, 1.0); CHKERRQ(ierr);
  ierr = IIMSetForceComponents(iim, MyForce ); CHKERRQ(ierr);

  IrregularNode *node;

  int i, k, len = 10;
  for( i = 0; i <= len; i++ ) {
    phi[0][0] = ( (len-i)*(-0.2) + i*(-5) ) / len;
    ierr = GridWrite(ls->phi, i); CHKERRQ(ierr);
    ierr = LevelSetUpdateIrregularNodeList(ls); CHKERRQ(ierr);
    ierr = IIMUpdateSurfaceQuantities( iim, ls); CHKERRQ(ierr);
    ierr = LevelSetWriteIrregularNodeList(ls, i); CHKERRQ(ierr);
    for (k = 0; k < ArrayLength(ls->irregularNodes); ++k) {
      ierr = ArrayGet(ls->irregularNodes,k,(void*)&node); CHKERRQ(ierr);
      if( node->axis == -1 ) continue;
      if( node->shift == CELL_CENTER ) {
        JumpVelocity( iim->mu, node, &jump, node->axis );
      } else {
        JumpVelocity( iim->mu, node, &jump, node->shift - U_FACE );
      }
      ierr = IIMLaplaceCorrection(iim, node, jump); CHKERRQ(ierr);
    }
    ierr = ArrayWrite(iim->val, i); CHKERRQ(ierr);
    ierr = ArrayWrite(iim->coor, i); CHKERRQ(ierr);
    ierr = ArraySetSize(iim->val,0); CHKERRQ(ierr);
    ierr = ArraySetSize(iim->coor,0); CHKERRQ(ierr);
  }


  ierr = LevelSetDestroy(ls); CHKERRQ(ierr);
  ierr = IIMDestroy(iim); CHKERRQ(ierr);
  ierr = DCellFinalize(); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

void MyForce( IrregularNode *n, void *ctx ) {
  n->F1 = 100*n->X.x;
}
