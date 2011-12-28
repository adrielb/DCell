#include "FluidField.h"

int main(int argc, char **args)
{
  PetscErrorCode  ierr;
  ierr = PetscInitialize(&argc, &args, (char *) 0, ""); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "Start\n"); CHKERRQ(ierr);
 
  int d1 = 32;
  int lsLen = 1;
  LevelSet lsets[lsLen];
  iCoor lsDim = {d1, d1, 0};
  ierr = LevelSetCreate(lsDim, &lsets[0]); CHKERRQ(ierr);
  ierr = LevelSetInitializeToStar(lsets[0]); CHKERRQ(ierr);
  
  ierr = PetscOptionsSetValue("-da_grid_x","32"); CHKERRQ(ierr);
  ierr = PetscOptionsSetValue("-da_grid_y","32"); CHKERRQ(ierr);
  FluidField f;
  ierr = FluidFieldCreate( &f ); CHKERRQ(ierr);
  
  PetscReal mu = 1;
  IIM iim;
  IIMCreate( &mu, lsDim, 32, &iim);
  IIMSetForceComponents(iim,ForceComponentNormalSurfaceTension,ForceComponentTangentialSurfaceTension);
  
  ierr = FluidFieldStep(f, iim, lsLen, lsets); CHKERRQ(ierr);
  /*
  VecWrite("p",f->p);
  VecWrite("px",f->px);
  VecWrite("py",f->py);
  VecWrite("u",f->u);
  VecWrite("v",f->v);
  */
  ierr = FluidFieldDestroy(f); CHKERRQ(ierr);
//  ierr = LevelSetDestroy(ls); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "End\n"); CHKERRQ(ierr);
  ierr = PetscFinalize(); CHKERRQ(ierr);
  return 0;
}