#include "FiberField.h"
/*#include "petscdmda.h"*/

#define napp 3

int main(int argc, char* args[])
{
  PetscErrorCode ierr;
  ierr = DCellInit(); CHKERRQ(ierr);
  
  PetscInt myapp[napp];
  PetscInt *mypetsc = 0;

  int rank;
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr); 
  if (rank == 0) {
    myapp[0] = 99;
    myapp[1] = 124;
    myapp[2] = 888;
  } 
  if (rank == 1 ) { 
    myapp[0] = 99;
    myapp[1] = 4;
    myapp[2] = 139000;
  }

  AO ao;
  /*ierr = AOCreateBasic(PETSC_COMM_WORLD, napp, myapp, mypetsc, &ao); CHKERRQ(ierr);*/
  ierr = AOCreateMapping(PETSC_COMM_WORLD, napp, myapp, mypetsc, &ao); CHKERRQ(ierr);
  ierr = AOView(ao, PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);

#define n 6
  PetscInt idx[n] = {4,99,888,139000,99,124};

  IS is;
  ierr = ISCreateGeneral(PETSC_COMM_WORLD, n, idx, PETSC_COPY_VALUES, &is); CHKERRQ(ierr);
  ierr = AOApplicationToPetscIS(ao, is); CHKERRQ(ierr);
  ierr = ISView(is, PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);

  ierr = ISDestroy(&is); CHKERRQ(ierr);
  ierr = AODestroy(&ao); CHKERRQ(ierr);
  ierr = DCellFinalize(); CHKERRQ(ierr);
  return 0;
}
