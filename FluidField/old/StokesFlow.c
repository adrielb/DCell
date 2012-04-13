#include "FluidField.h"
#include "ImmersedInterfaceMethod.h"

int main(int argc, char **args)
{
  PetscErrorCode  ierr;
  
  ierr = PetscInitialize(&argc, &args, (char *) 0, ""); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "Start\n"); CHKERRQ(ierr);
  
  ierr = PetscOptionsSetValue("-da_grid_x","63"); CHKERRQ(ierr);
  ierr = PetscOptionsSetValue("-da_grid_y","63"); CHKERRQ(ierr);

  IIM iim;
  IIMCreate(12, &iim);
  IIMSetForceComponents(iim,ForceComponentNormalSimple,ForceComponentTangentialSimple);
  
  PetscInt d1=63,d2=d1;
  LevelSet2D ls;
  CreateLevelSet2D(d1,d2,&ls);
  LevelSetInitializeToCircle(ls, PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE);
  IIMUpdateSurfaceQuantities2D( iim, ls );
  IIMUpdateSurfaceDerivatives2D( iim, ls );
  IrregularNodeListWrite(ls->irregularNodes,0);
  
  FluidField f;
  ierr = FluidFieldCreate(&f); CHKERRQ(ierr);
  DALocalInfo info;
  ierr = DAGetLocalInfo(f->da,&info); CHKERRQ(ierr);
  PetscReal dx = f->d.x, dy = f->d.y;

  PetscReal **px, **py;
  ierr = DAVecGetArray(f->da, f->px, &px); CHKERRQ(ierr);
  ierr = DAVecGetArray(f->da, f->py, &py); CHKERRQ(ierr);
  //[ px, py ] = grad p
  ierr = Gradient2D(info,dx,dy,f->p, f->px, f->py); CHKERRQ(ierr);
    
  //[ px, py] += pc
  ierr = IIMSimplePressureGradientCorrection2D(iim, ls, f->d, px, py); CHKERRQ(ierr);
  ierr = VecWrite("px",f->px); CHKERRQ(ierr);
  ierr = VecWrite("py",f->py); CHKERRQ(ierr);
  
  // px /= mu
  ierr = VecScale(f->px, 1/f->mu); CHKERRQ(ierr);
// px += uc
  ierr = IIMSimpleCorrection2D(iim, ls, JumpConditionVelocity2D_X, f->d, px); CHKERRQ(ierr);

//L(u) = ( px + pc + uc ) / mu
  ierr = KSPSolve(f->kspU,f->px,f->u); CHKERRQ(ierr);
  ierr = VecWrite("u",f->u); CHKERRQ(ierr);
  
// py += vc
  ierr = IIMSimpleCorrection2D(iim, ls, JumpConditionVelocity2D_Y, f->d, py); CHKERRQ(ierr);
// py /= mu
  ierr = VecScale(f->py, 1/f->mu); CHKERRQ(ierr);
//L(u) = ( py + pc + vc ) / mu
  ierr = KSPSolve(f->kspV,f->py,f->v); CHKERRQ(ierr);
  ierr = VecWrite("v",f->v); CHKERRQ(ierr);

//div = ux + vy
  ierr = Divergence(info,dx,dy,f->u,f->v, f->div); CHKERRQ(ierr);
  ierr = VecWrite("divV",f->div); CHKERRQ(ierr); 

  VecZeroEntries(f->div);
  PetscReal **div;
  ierr = DAVecGetArray(f->da, f->div, &div); CHKERRQ(ierr);
  ierr = IIMSimpleCorrection2D(iim, ls, JumpConditionPressure2D, f->d, div); CHKERRQ(ierr);
// Projection step?
  ierr = KSPSolve(f->kspPhi,f->div,f->p); CHKERRQ(ierr);
  ierr = VecWrite("p",f->p); CHKERRQ(ierr);
  
  ierr = Gradient2D(info,dx,dy,f->p,f->px,f->py); CHKERRQ(ierr);
  //[ px, py] += pc
  ierr = IIMSimplePressureGradientCorrection2D(iim, ls, f->d, px, py); CHKERRQ(ierr);
  ierr = VecAXPY(f->u,-1,f->px); CHKERRQ(ierr);
  ierr = VecAXPY(f->v,-1,f->py); CHKERRQ(ierr);
  ierr = VecWrite("uc",f->u); CHKERRQ(ierr);
  ierr = VecWrite("vc",f->v); CHKERRQ(ierr);
  
  ierr = Divergence(info,dx,dy,f->u,f->v, f->div); CHKERRQ(ierr);
  ierr = VecWrite("divC",f->div); CHKERRQ(ierr);

  ierr = FluidFieldDestroy(f); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "End\n"); CHKERRQ(ierr);
  ierr = PetscFinalize(); CHKERRQ(ierr);
  return 0;
}
