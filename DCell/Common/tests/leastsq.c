#include "Common.h"

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc, char **args)
{
  int i;
  PetscErrorCode ierr;
  ierr = DCellInit(); CHKERRQ(ierr);
  int Np = 10;
  PetscBool is2D = PETSC_TRUE;
  LeastSq lsq;
  ierr = LeastSqCreate( Np, is2D, &lsq); CHKERRQ(ierr);
  PetscReal *s, *g;
  ierr = LeastSqGetVecs(lsq,&s,0,&g,0); CHKERRQ(ierr);
  ierr = LeastSqSetNumPoints(lsq,Np); CHKERRQ(ierr);
  for (i = 0; i < Np; ++i) {
    s[i] = i - Np/2;
    g[i] = s[i]*s[i];
  }
  ierr = LeastSqSolve(lsq); CHKERRQ(ierr);
  for (i = 0; i < Np; i++) {
    printf("g[%d] = %f\n", i, g[i] );
  }

  ierr = LeastSqDestroy(lsq); CHKERRQ(ierr);
	ierr = DCellFinalize(); CHKERRQ(ierr);
  return 0;
}
