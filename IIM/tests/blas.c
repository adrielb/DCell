#include "ImmersedInterfaceMethod.h"

void drot_(int *len, double* x, int *incx, double* y, int *incy, double *c, double *s);

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc, char **args)
{
  PetscErrorCode  ierr;
  PetscFunctionBegin;
  ierr = PetscInitialize(&argc, &args, (char *) 0, ""); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "Start %s\n", __FILE__); CHKERRQ(ierr);

  int len = 3, incx = 1, incy = 1;
  double x[] = { 1,2,3 }, c = .5;
  double y[] = { 1,2,3 }, s = .5;


  drot_(&len,x,&incx,y,&incy,&c,&s);

  for( int i = 0; i < len; i++ ){
    printf("{%f, %f},",x[i],y[i]);
  }
  printf("\n");

  ierr = PetscPrintf(PETSC_COMM_WORLD, "End %s\n", __FILE__); CHKERRQ(ierr);
  ierr = PetscFinalize(); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
