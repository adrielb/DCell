#include "ImmersedInterfaceMethod.h"

void InterfacialForceSpin( IrregularNode *n, void *ctx );

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc, char **args)
{
  PetscErrorCode  ierr;
  PetscFunctionBegin;
  ierr = DCellInit(); CHKERRQ(ierr);

  PetscReal dx = 2;
  Coor dh = {dx, dx, 0};
  Coor center = {12,12,0};
  PetscReal radius = 9;

  LevelSet ls;
  ierr = LevelSetInitializeToCircle(dh,center, radius, &ls); CHKERRQ(ierr);
//  ierr = LevelSetInitializeToStar2D(dh,center,radius,0.5*radius, 5,&ls); CHKERRQ(ierr);
  ierr = GridWrite(ls->phi,0); CHKERRQ(ierr);
  ierr = LevelSetWriteIrregularNodeList(ls,-1); CHKERRQ(ierr);

  PetscReal mu=1;
  IIM iim;
  ierr = IIMCreate(ls->phi->is2D,&mu,1,20,dh,&iim); CHKERRQ(ierr);
  ierr = IIMSetForceComponents(iim,InterfacialForceSpin); CHKERRQ(ierr);
  ierr = IIMUpdateSurfaceQuantities(iim,ls); CHKERRQ(ierr);
  ierr = LevelSetWriteIrregularNodeList(ls,1); CHKERRQ(ierr);

  DM da;
  int dof=3;
  iCoor n = ls->phi->n;
  n.x *= 0.5;
  n.y *= 0.5;
  printf("%d x %d\n", n.x, n.y );
  ierr = DMDACreate2d(PETSC_COMM_WORLD,DMDA_BOUNDARY_NONE,DMDA_BOUNDARY_NONE,DMDA_STENCIL_BOX,
            n.x,n.y,PETSC_DECIDE,PETSC_DECIDE,dof,1,0,0,&da); CHKERRQ(ierr);

  Vec vecRHS;
  ierr = DMCreateGlobalVector(da,&vecRHS); CHKERRQ(ierr);

  int ga;
  ierr = GACreate(da,&ga); CHKERRQ(ierr);

  ierr = IIMUpdateRHS( iim, ls, ga ); CHKERRQ(ierr);

  ierr = GAGetVec(ga,vecRHS); CHKERRQ(ierr);

  ierr = VecWrite(vecRHS,"rhs",0); CHKERRQ(ierr);

  ierr = DMDestroy(&da); CHKERRQ(ierr);
  ierr = VecDestroy(&vecRHS); CHKERRQ(ierr);
  ierr = LevelSetDestroy(ls); CHKERRQ(ierr);
  ierr = IIMDestroy(iim); CHKERRQ(ierr);
  ierr = DCellFinalize(); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

void InterfacialForceSpin( IrregularNode *n, void *ctx )
{
  n->f1 = 0;
  n->f2 = 1;
}
