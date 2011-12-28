#include "FluidField.h"


#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc, char **args)
{
  PetscErrorCode  ierr;
  ierr = PetscInitialize(&argc, &args, (char *) 0, ""); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "Start %s\n", __FILE__); CHKERRQ(ierr);

  FluidField f;
  ierr = FluidFieldCreate( PETSC_COMM_WORLD, &f); CHKERRQ(ierr);
  ierr = FluidFieldWrite(f,0); CHKERRQ(ierr);
  ierr = FluidFieldDestroy(f); CHKERRQ(ierr);

  ierr = PetscPrintf(PETSC_COMM_WORLD, "End %s\n", __FILE__); CHKERRQ(ierr);
  ierr = PetscFinalize(); CHKERRQ(ierr);
}
