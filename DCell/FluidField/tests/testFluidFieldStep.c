#include "FluidField.h"

PetscErrorCode RHSFunc( PetscReal dx, int t, FluidField f );

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc, char **args)
{
  PetscErrorCode  ierr;
  ierr = DCellInit(); CHKERRQ(ierr);

  PetscReal d1 = 256, dx = 1./d1;
  iCoor size = {d1,d1,0};

  FluidField fluid;
  ierr = FluidFieldCreate(PETSC_COMM_WORLD, &fluid); CHKERRQ(ierr);
  ierr = FluidFieldSetDims(fluid,size); CHKERRQ(ierr);
  ierr = FluidFieldSetDx(fluid,dx); CHKERRQ(ierr);
  ierr = FluidFieldSetup(fluid); CHKERRQ(ierr);

  for (int t = 0; t < 4; ++t) {
    ierr = VecSet( fluid->rhs, 0); CHKERRQ(ierr);
    ierr = RHSFunc( dx, t, fluid ); CHKERRQ(ierr);
    ierr = FluidFieldSolve(fluid); CHKERRQ(ierr);
//    ierr = VecWrite(fluid->vel,"vel",t); CHKERRQ(ierr);
  }

  ierr = FluidFieldDestroy(fluid); CHKERRQ(ierr);
  ierr = DCellFinalize(); CHKERRQ(ierr);
}

PetscErrorCode RHSFunc( PetscReal dx, int t, FluidField f )
{
  int xs, xm, ys, ym;
  PetscReal ***rhs;
  PetscReal s = t / 100.;
  PetscErrorCode ierr;

  ierr = DAGetCorners(f->daV,&xs,&ys,0,&xm,&ym,0); CHKERRQ(ierr);
  ierr = DAVecGetArrayDOF(f->daV,f->rhs,&rhs); CHKERRQ(ierr);
  for (int j = ys; j < ys+ym; ++j) {
    for (int i = xs; i < xs+xm; ++i) {
      if ( 0.25+s <= i*dx && i*dx <= 0.75+s &&
           0.25+s <= j*dx && j*dx <= 0.75+s ) {
        rhs[j][i][0] = 0.001 / (dx*dx);
        rhs[j][i][1] = 0.001 / (dx*dx);
      }
    }
  }
  ierr = DAVecRestoreArrayDOF(f->daV,f->rhs,&rhs); CHKERRQ(ierr);

  return 0;
}
