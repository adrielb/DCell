#include "FiberField.h"

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc, char **args) {
  PetscErrorCode ierr;
  ierr = DCellInit(); CHKERRQ(ierr);
  MPI_Comm comm = PETSC_COMM_WORLD;
  
  int tag = 0; 
  MPI_Request reqSend;
  MPI_Request reqRecv;

  int size;
  MPI_Comm_size( comm, &size);
  ierr = PetscPrintf( comm, "Comm size: %d\n", size); CHKERRQ(ierr);
  int rank; 
  MPI_Comm_rank( comm, &rank ); 

  int data = rank + 10;
  ierr = MPI_Isend( &data, 1, MPI_INT, rank, tag, comm, &reqSend ); CHKERRQ(ierr);

  int buf;
  ierr = MPI_Irecv( &buf, 1, MPI_INT, rank, tag, comm, &reqRecv ); CHKERRQ(ierr);

  ierr = PetscSynchronizedPrintf( comm, "[%d] %d\n",rank, buf); CHKERRQ(ierr);
  ierr = PetscSynchronizedFlush( comm ); CHKERRQ(ierr);

  ierr = DCellFinalize(); CHKERRQ(ierr);
  return 0;
}
