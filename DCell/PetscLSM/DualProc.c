#include "petscksp.h"

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc, char **args)
{
  PetscErrorCode  ierr;
  PetscInitialize(&argc, &args, (char *)0, "");
  PetscPrintf(PETSC_COMM_WORLD, "Start\n");
  
  PetscInt rank;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  PetscInt size;
  MPI_Comm_size(PETSC_COMM_WORLD, &size);
  MPI_Comm commSMP;
  MPI_Comm_split(PETSC_COMM_WORLD, rank/2, rank, &commSMP);
  char name[MPI_MAX_PROCESSOR_NAME];
  int namelen;
  MPI_Get_processor_name(name, &namelen);
  
  PetscPrintf(commSMP, "%s, %d - %d, %d\n", name, rank, size, rank/2);
  
  ierr = PetscFinalize(); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}