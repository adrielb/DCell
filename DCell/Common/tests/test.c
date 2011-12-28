#include "petsc.h"

static inline PetscReal diff( PetscReal a1, PetscReal a2 ){
  return a1 - a2;
}

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc, char **args)
{
  PetscErrorCode  ierr;
  PetscFunctionBegin;
  ierr = PetscInitialize(&argc, &args, (char *) 0, ""); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "Start %s\n", __FILE__); CHKERRQ(ierr);
  
  int n = 5e5;
  PetscReal *a1,*a2;
  ierr = PetscMalloc(sizeof(PetscReal)*n,&a1); CHKERRQ(ierr);
  ierr = PetscMalloc(sizeof(PetscReal)*n,&a2); CHKERRQ(ierr);

  int j;
  printf("Without Function Call\n");
  for (j = 0; j < 10; ++j) {
    int i;
    PetscLogDouble t1, t2;
    ierr = PetscGetTime(&t1); CHKERRQ(ierr);
    for( i = 0; i < n; i++ )
    {
      if( (a1[i] - a2[i]) < 0 ) {
        a1[i] = i;
      }
      a2[i] = i;
    }
    ierr = PetscGetTime(&t2); CHKERRQ(ierr);
    printf("TIME: %f\n", t2-t1);
  }

  printf("With Function Call\n");
  for (j = 0; j < 10; ++j) {
    int i;
    PetscLogDouble t1, t2;
    ierr = PetscGetTime(&t1); CHKERRQ(ierr);
    for( i = 0; i < n; i++ )
    {
      if( diff(a1[i], a2[i]) < 0 ) {
        a1[i] = i;
      }
      a2[i] = i;
    }
    ierr = PetscGetTime(&t2); CHKERRQ(ierr);
    printf("TIME: %f\n", t2-t1);
  }

  for( j = 0; j < 10; j++ )
  {
    printf("j: %d\n", j);
    if( j==9 ) break;
  }
  printf("J: %d\n", j);

  ierr = PetscPrintf(PETSC_COMM_WORLD, "End %s\n", __FILE__); CHKERRQ(ierr);
  ierr = PetscFinalize(); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
