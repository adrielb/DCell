#include "FluidField.h"

void DetermineTimeStep(FluidField f, PetscReal *dt);

int main(int argc, char **args)
{
  PetscErrorCode  ierr;
  ierr = PetscInitialize(&argc, &args, (char *) 0, ""); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "Start\n"); CHKERRQ(ierr);
  
                                   int d1 = 9;
  ierr = PetscOptionsSetValue("-da_grid_x","9"); CHKERRQ(ierr);
  ierr = PetscOptionsSetValue("-da_grid_y","9"); CHKERRQ(ierr);
//  ierr = PetscOptionsSetValue("-da_grid_z","9"); CHKERRQ(ierr);
  
  FluidField f;
  ierr = FluidFieldCreate( 1, 1, &f); CHKERRQ(ierr);
  
  ierr = MatWrite("U",f->matU); CHKERRQ(ierr);
  ierr = MatWrite("V",f->matV); CHKERRQ(ierr);
  if( f->matW ) ierr = MatWrite("W",f->matW); CHKERRQ(ierr);
  ierr = MatWrite("P",f->matP); CHKERRQ(ierr);
  
  ierr = FluidFieldDestroy(f); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "End\n"); CHKERRQ(ierr);
  ierr = PetscFinalize(); CHKERRQ(ierr);
  return 0;
}