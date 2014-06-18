#include "FiberField.h"

#undef __FUNCT__
#define __FUNCT__ "FiberFieldSetFluidVelocityEvaluator"
PetscErrorCode FiberFieldSetFluidVelocityEvaluator( FiberField f, void *func )
{
  PetscFunctionBegin;

  f->EvaluateFluidVelocity = func;

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "FiberFieldEvaluateFluidVelocity"
PetscErrorCode FiberFieldEvaluateFluidVelocity( FiberField f, PetscReal t, Vec x )
{
  PetscErrorCode ierr;

  PetscFunctionBegin;

  ierr = f->EvaluateFluidVelocity( f, t, x); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "FiberField_ZeroFluidVelocity"
PetscErrorCode FiberField_ZeroFluidVelocity( FiberField f, PetscReal t, Vec x )
{
  PetscErrorCode ierr;

  PetscFunctionBegin;

  ierr = VecZeroEntries( f->vf ); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "FiberField_CircularFluidVelocity"
PetscErrorCode FiberField_CircularFluidVelocity( FiberField f, PetscReal t, Vec X )
{
  int i;
  int len;
  PetscReal *x;
  PetscReal *vf;
  PetscErrorCode ierr;

  PetscFunctionBegin;

  ierr = VecGetLocalSize( X, &len); CHKERRQ(ierr);
  ierr = VecGetArray( X, &x); CHKERRQ(ierr);
  ierr = VecGetArray( f->vf, &vf); CHKERRQ(ierr);
  for (i = 0; i < len; i+=3) {
    vf[i+0] = -x[i+2]; // V.x = -X.z
    vf[i+1] = 0;       // V.y =  0 
    vf[i+2] =  x[i+0]; // V.z =  X.x
  }
  ierr = VecRestoreArray( X, &x); CHKERRQ(ierr);
  ierr = VecRestoreArray( f->vf, &vf); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}
