#include "petsc.h"
#include "ImmersedInterfaceMethod.h"

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc, char **args)
{
  PetscErrorCode  ierr;
  PetscFunctionBegin;
  ierr = PetscInitialize(&argc, &args, (char *) 0, ""); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "Start\n"); CHKERRQ(ierr);
  
  IIM iim;
  IIMCreate(12, &iim);
  IIMSetForceComponents(iim,ForceComponentNormalSimple,ForceComponentTangentialSimple);
  
  iCoor s = {63, 63, 0};
  LevelSet ls;
  LevelSetCreate(s,&ls);
  LevelSetInitializeToCircle(ls, PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE);
  
  IrregularNodeListWrite(ls->irregularNodes, 0);
  
  Coor dh = {1,1,1};
  Grid px, py, pz;
  GridCreate(s, &px);
  GridCreate(s, &py);
  IIMUpdateSurfaceQuantities2D( iim, ls );
  IIMUpdateSurfaceDerivatives2D( iim, ls);
  IIMSimplePressureGradientCorrection2D(iim, ls, dh, px->v2,py->v2);
  WriteVector("px",px->v);
  WriteVector("py",py->v);
  
  Grid rhs;
  GridCreate(s, &rhs);
  IIMSimpleCorrection2D(iim, ls, JumpConditionPressure2D, dh, rhs->v2);
  WriteVector("pc",rhs->v);
  VecZeroEntries(rhs->v);
  
  IrregularNodeListUpdate(-1,0,ls); 
  IIMUpdateSurfaceQuantities2D(iim,ls);
  IIMUpdateSurfaceDerivatives2D(iim,ls);
  IIMSimpleCorrection2D(iim, ls, JumpConditionVelocity2D_X, dh, rhs->v2);
  WriteVector("uc",rhs->v);
  VecZeroEntries(rhs->v);
  
  IrregularNodeListUpdate(0,-1,ls);
  IIMUpdateSurfaceQuantities2D(iim,ls);
  IIMUpdateSurfaceDerivatives2D(iim,ls); 
  IIMSimpleCorrection2D(iim, ls, JumpConditionVelocity2D_Y, dh, rhs->v2);
  WriteVector("vc",rhs->v);
  
  LevelSetDestroy(ls);
  IIMDestroy(iim);
  
  ierr = PetscPrintf(PETSC_COMM_WORLD, "End\n"); CHKERRQ(ierr);
  ierr = PetscFinalize(); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}