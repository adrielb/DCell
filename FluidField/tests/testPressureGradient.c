#include "FluidField.h"

int main(int argc, char **args)
{
  PetscErrorCode  ierr;
  ierr = PetscInitialize(&argc, &args, (char *) 0, ""); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "Start\n"); CHKERRQ(ierr);
  
  int d1=63, d2=d1;
  ierr = PetscOptionsSetValue("-da_grid_x","63"); CHKERRQ(ierr);
  ierr = PetscOptionsSetValue("-da_grid_y","63"); CHKERRQ(ierr);
    
  LevelSet2D ls;
  CreateLevelSet2D(d1,d2,&ls);
  LevelSetInitializeToCircle(ls, PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE);
  
  IIM iim;
  IIMCreate(12, &iim);
  IIMSetForceComponents(iim,ForceComponentNormalSimple,ForceComponentTangentialSimple);
  
  FluidField f;
  FluidFieldCreate(&f);
    
  Grid2D px, py;
  CreateGrid2D(d1,d2,&px);
  CreateGrid2D(d1,d2,&py);
  IIMUpdateSurfaceQuantities2D( iim, ls );
  IIMUpdateSurfaceDerivatives2D( iim, ls);
  IIMSimplePressureGradientCorrection2D(iim,ls, f->d, px->v2,py->v2);
  VecWrite("pxc",px->v);
  VecWrite("pyc",py->v);
  
  VecZeroEntries(f->div);
  PetscReal **div;
  ierr = DAVecGetArray(f->da, f->div, &div); CHKERRQ(ierr);
  ierr = IIMSimpleCorrection2D(iim, ls, JumpConditionPressure2D, f->d, div); CHKERRQ(ierr);
  ierr = KSPSolve(f->kspPhi,f->div,f->p); CHKERRQ(ierr);
  VecWrite("p",f->p);
  DALocalInfo info;
  ierr = DAGetLocalInfo(f->da,&info); CHKERRQ(ierr);
  ierr = Gradient2D(info,f->d.x, f->d.y, f->p, f->px, f->py); CHKERRQ(ierr);
  VecWrite("px",f->px);
  VecWrite("py",f->py);
  
  FluidFieldDestroy(f);
  DestroyLevelSet2D(ls);
  IIMDestroy(iim);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "End\n"); CHKERRQ(ierr);
  ierr = PetscFinalize(); CHKERRQ(ierr);
  return 0;
}