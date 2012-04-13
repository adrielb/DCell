#include "LevelSetMethod.h"

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc, char **args)
{
  PetscErrorCode  ierr;
  PetscFunctionBegin;
  ierr = PetscInitialize(&argc, &args, (char *) 0, ""); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "Start %s\n", __FILE__); CHKERRQ(ierr);

  int d1 = 512;
  PetscReal dx = 1./(d1-1);
  Coor dh = {dx, dx, dx};
  iCoor pos = {0,0,0};
  iCoor size = {d1,d1,d1};
  Grid phi;
  ierr = GridCreate(dh,pos,size,1,&phi); CHKERRQ(ierr);
  PetscReal **phi2D, ***phi3D;
  ierr = GridGet(phi,&phi2D); CHKERRQ(ierr);
  ierr = GridGet(phi,&phi3D); CHKERRQ(ierr);

  int i,j,k;
  PetscLogDouble t1,t2;

  for (j = 0; j < 6; ++j) {
    for (i = 0; i < 3; ++i) {
      printf("%d", STAR[j][i]);
    }
    printf("\n");
  }

  ierr = PetscGetTime(&t1); CHKERRQ(ierr);
  for (k = 0; k < d1; ++k) {
    for (j = 0; j < d1; ++j) {
      for (i = 0; i < d1; ++i) {
        phi3D[k][j][i] = i*j*k;
      }
    }
  }
  ierr = PetscGetTime(&t2); CHKERRQ(ierr);
  printf("time: %f\n", t2-t1);

  ierr = PetscGetTime(&t1); CHKERRQ(ierr);
  for (k = 0; k < (phi->is2D ? 1 : d1); ++k) {
    for (j = 0; j < d1; ++j) {
      for (i = 0; i < d1; ++i) {
        if( phi->is2D ) {
          phi2D[j][i] = i*j;
        } else {
          phi3D[k][j][i] = i*j*k;
        }
      }
    }
  }
  ierr = PetscGetTime(&t2); CHKERRQ(ierr);
  printf("time: %f\n", t2-t1);

  const PetscTruth is2D = phi->is2D;
  ierr = PetscGetTime(&t1); CHKERRQ(ierr);
  for (k = 0; k < (is2D ? 1 : d1); ++k) {
    for (j = 0; j < d1; ++j) {
      for (i = 0; i < d1; ++i) {
        if( is2D ) {
          phi2D[j][i] = i*j;
        } else {
          phi3D[k][j][i] = i*j*k;
        }
      }
    }
  }
  ierr = PetscGetTime(&t2); CHKERRQ(ierr);
  printf("time: %f\n", t2-t1);

  ierr = GridDestroy(phi); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "End %s\n", __FILE__); CHKERRQ(ierr);
  ierr = PetscFinalize(); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
