#include "DWorld.h"

static PetscReal mu = 0.02;

void MyInterfacialForce(IrregularNode *n, void *context )
{
  n->F1 = 0;
  n->F2 = 10*mu;
//  n->F1 = n->k;
//  n->F2 = 0;
}

#undef __FUNCT__
#define __FUNCT__ "UpdateFluidFieldRHS"
PetscErrorCode UpdateFluidFieldRHS( DCell dcell, IIM iim, int ga, PetscReal t ) {
  PetscErrorCode ierr;
  PetscFunctionBegin
  ierr = IIMSetForceComponents(iim, InterfacialForceSurfaceTension ); CHKERRQ(ierr);
  ierr = IIMSetForceContext(iim, dcell); CHKERRQ(ierr);
  ierr = IIMUpdateRHS(iim, dcell->lsPlasmaMembrane, ga); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc, char **args) {
  PetscErrorCode ierr;
  ierr = DCellInit(); CHKERRQ(ierr);

  PetscReal dx = 1./128.;
  Coor len = {2,2,0};
  Coor dh  = {dx,dx,0};
  iCoor size = {len.x/dx,len.y/dx,0};
  printf("MX = %d;\n", size.x);
  printf("MY = %d;\n", size.y);

  FluidField fluid;
  ierr = FluidFieldCreate(PETSC_COMM_WORLD, &fluid);  CHKERRQ(ierr);
  ierr = FluidFieldSetDims(fluid,size); CHKERRQ(ierr);
  ierr = FluidFieldSetDx(fluid,dx); CHKERRQ(ierr);
  ierr = FluidFieldSetViscosity(fluid,mu); CHKERRQ(ierr);
  ierr = FluidFieldSetup(fluid); CHKERRQ(ierr);

  IIM iim;
  ierr = IIMCreate(!fluid->is3D,&fluid->mu,1.0,32,dh,&iim); CHKERRQ(ierr);

  Coor center = (Coor){1,1,0};
  PetscReal radius = 0.5001;
  LevelSet ls;
//  ierr = LevelSetInitializeToStar2D(fluid->dh,center,radius,0.5*radius, 5,&ls); CHKERRQ(ierr);
  ierr = LevelSetInitializeToCircle(dh, center, radius, &ls); CHKERRQ(ierr);
  ierr = LevelSetWriteIrregularNodeList( ls->irregularNodes, 0 ); CHKERRQ(ierr);
  ierr = GridWrite(ls->phi,0); CHKERRQ(ierr);

  ierr = IIMSetForceComponents(iim, MyInterfacialForce); CHKERRQ(ierr);
  ierr = IIMUpdateRHS(iim, ls, fluid->ga); CHKERRQ(ierr);
  ierr = GAGetVec(fluid->ga,fluid->rhs); CHKERRQ(ierr);

  ierr = FluidFieldSolve( fluid ); CHKERRQ(ierr);
  printf("AGAIN!!!\n");
  ierr = FluidFieldSolve( fluid ); CHKERRQ(ierr);
  ierr = FluidFieldWrite( fluid, 0); CHKERRQ(ierr);
  ierr = MatWrite(fluid->mat,"mat",0); CHKERRQ(ierr);

  ierr = DCellFinalize(); CHKERRQ(ierr);
  return 0;
}
