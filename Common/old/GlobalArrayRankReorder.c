#include "petsc.h"

PetscErrorCode ReorderRanks();

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc, char **args)
{
  PetscErrorCode  ierr;
  PetscFunctionBegin;
//  ierr = ReorderRanks();
  ierr = PetscInitialize(&argc, &args, (char *) 0, ""); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "Start %s\n", __FILE__); CHKERRQ(ierr);

  int xprocs, yprocs, zprocs;
  PetscTruth flgFound;
  ierr = PetscOptionsGetInt(PETSC_NULL,"-da_processors_x",&xprocs,&flgFound); CHKERRQ(ierr);
  ierr = PetscOptionsGetInt(PETSC_NULL,"-da_processors_y",&yprocs,&flgFound); CHKERRQ(ierr);
  ierr = PetscOptionsGetInt(PETSC_NULL,"-da_processors_z",&zprocs,&flgFound); CHKERRQ(ierr);

//  ierr = PetscSynchronizedPrintf(PETSC_COMM_WORLD,"xprocs: %d\n", xprocs); CHKERRQ(ierr);
//  ierr = PetscSynchronizedFlush(PETSC_COMM_WORLD); CHKERRQ(ierr);



  ierr = PetscPrintf(PETSC_COMM_WORLD, "End %s\n", __FILE__); CHKERRQ(ierr);
  ierr = PetscFinalize(); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
/*
PetscErrorCode ReorderRanks()
{
  MPI_Comm NewComm;
  int MPI_Rank, NewRank, x,y,z;

  // get rank from MPI ordering:
  MPI_Comm_rank(MPI_COMM_WORLD, &MPI_Rank);

  // calculate coordinates of cpus in MPI ordering:
  x = MPI_rank / (z_procs*y_procs);
  y = (MPI_rank % (z_procs*y_procs)) / z_procs;
  z = (MPI_rank % (z_procs*y_procs)) % z_procs;

  // set new rank according to PETSc ordering:
  NewRank = z*y_procs*x_procs + y*x_procs + x;

  // create communicator with new ranks according to PETSc ordering:
  MPI_Comm_split(PETSC_COMM_WORLD, 1, NewRank, &NewComm);

  // override the default communicator (was MPI_COMM_WORLD as default)
  PETSC_COMM_WORLD = NewComm;
}
*/
