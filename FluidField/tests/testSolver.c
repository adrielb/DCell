#include "DWorld.h"

typedef struct _MyCell {
  struct _DCell dcell;
  PetscReal K;  // magnitude of surface tension
} *MyCell;

#undef __FUNCT__
#define __FUNCT__ "MyCellCreate"
PetscErrorCode MyCellCreate( LevelSet ls, MyCell *mycell )
{
  MyCell cell;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscNew( struct _MyCell, &cell); CHKERRQ(ierr);
  ierr = DCellSetup( ls, (DCell)cell ); CHKERRQ(ierr);

  *mycell = cell;
  PetscFunctionReturn(0);
}

void InterfacialForceAdhesionFeedback(IrregularNode *n, void *context )
{
  const MyCell c = (MyCell)context;

  n->F1 = -c->K * n->k;
  n->F2 =  0;
}

PetscErrorCode UpdateFluidFieldRHS( DCell dcell, IIM iim, int ga, PetscReal t ) {
  PetscErrorCode ierr;
  PetscFunctionBegin
  ierr = IIMSetForceComponents(iim, InterfacialForceAdhesionFeedback ); CHKERRQ(ierr);
  ierr = IIMSetForceContext(iim, dcell); CHKERRQ(ierr);
  ierr = IIMUpdateRHS(iim, dcell->lsPlasmaMembrane, ga); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

int main(int argc, char **args) {
  PetscErrorCode ierr;
  ierr = DCellInit(); CHKERRQ(ierr);

  int d1 = 32;
  PetscReal dx = 1./d1;
  FluidField fluid;
  ierr = FluidFieldCreate(PETSC_COMM_WORLD, &fluid);  CHKERRQ(ierr);
  ierr = FluidFieldSetDims(fluid,(iCoor){2*d1,d1}); CHKERRQ(ierr);
  ierr = FluidFieldSetDx(fluid,dx); CHKERRQ(ierr);
  ierr = FluidFieldSetup(fluid); CHKERRQ(ierr);

  DWorld world;
  ierr = DWorldCreate(fluid, &world); CHKERRQ(ierr);

  Coor center = (Coor){1,.4,0};
  PetscReal radius = 0.4-dx/2.;
  LevelSet ls;
  ierr = LevelSetInitializeToCircle(fluid->dh,center,radius,&ls); CHKERRQ(ierr);
  MyCell cell;
  ierr = MyCellCreate( ls, &cell ); CHKERRQ(ierr);

  cell->K = 100;
  cell->dcell.UpdateFluidFieldRHS = UpdateFluidFieldRHS;

  ierr = DWorldAddDCell( world, cell ); CHKERRQ(ierr);
  world->timax = 5000;
  world->dtmax = 1e4;
  world->CFL = 0.1;
  world->tend = 1e3;

  ierr = DWorldSimulate(world); CHKERRQ(ierr);
  ierr = DWorldDestroy(world); CHKERRQ(ierr);
  ierr = DCellFinalize(); CHKERRQ(ierr);
  return 0;
}
