#include "DWorld.h"

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc, char **args) {

  PetscErrorCode ierr;
  ierr = DCellInit(); CHKERRQ(ierr);

  int d1 = 32;
  FluidField fluid;
  ierr = FluidFieldCreate(PETSC_COMM_WORLD, &fluid);  CHKERRQ(ierr);
  ierr = FluidFieldSetDims(fluid,(iCoor){2*d1,d1}); CHKERRQ(ierr);
  ierr = FluidFieldSetDx(fluid,1./d1); CHKERRQ(ierr);
  ierr = FluidFieldSetup(fluid); CHKERRQ(ierr);

  DWorld world;
  ierr = DWorldCreate(fluid, &world); CHKERRQ(ierr);

  LevelSet ls;
  ierr = LevelSetInitializeToCircle((Coor){1,1,1},(Coor){0,0,0},5,&ls); CHKERRQ(ierr);

  DCell dcell;
  ierr = DCellCreate(ls, &dcell); CHKERRQ(ierr);
  ierr = DWorldAddDCell(world, dcell ); CHKERRQ(ierr);
  ierr = DWorldDestroy(world); CHKERRQ(ierr);

	ierr = DCellFinalize(); CHKERRQ(ierr);
  return 0;
}
