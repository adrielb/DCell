#include "petsc.h"
#include <stdio.h>

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc, char **args)
{
  PetscErrorCode  ierr;
  ierr = PetscInitialize(&argc, &args, (char *) 0, ""); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "Start %s\n", __FILE__); CHKERRQ(ierr);
  MPI_Comm comm = PETSC_COMM_WORLD;

  int rank;
  MPI_Comm_rank(comm,&rank);
  char hostname[256];
  gethostname(hostname,255);
  printf("%d: %s\n",rank,hostname);
  PetscBarrier(0);

  ierr = PetscPrintf(PETSC_COMM_WORLD, "End %s\n", __FILE__); CHKERRQ(ierr);
  ierr = PetscFinalize(); CHKERRQ(ierr);
}
