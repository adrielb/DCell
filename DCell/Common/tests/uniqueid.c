#include "Common.h"

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc, char **args) {
  PetscErrorCode ierr;
  ierr = DCellInit(); CHKERRQ(ierr);
  MPI_Comm comm = MPI_COMM_WORLD;
  int size, rank, i;
  MPI_Comm_size(comm,&size);
  MPI_Comm_rank(comm,&rank);

  ierr = PetscPrintf(comm,"\nCreating\n"); CHKERRQ(ierr);
  UniqueID uid;
  ierr = UniqueIDCreate(&uid); CHKERRQ(ierr);



  UniqueIDType id;
  ierr = PetscPrintf(comm,"\nLocal Count == 0\n"); CHKERRQ(ierr);
  for (i = 0; i < 2*size; ++i) {
    ierr = UniqueIDGenerate(uid,&id); CHKERRQ(ierr);
    ierr = PetscSynchronizedPrintf(comm,"%d: %d\n", rank, id); CHKERRQ(ierr);
  }
  ierr = PetscSynchronizedFlush(comm); CHKERRQ(ierr);

  int maxID = 5*size*size;
  ierr = PetscPrintf(comm,"\nSetting max id = %d\n", maxID); CHKERRQ(ierr);
  ierr = UniqueIDSetStartCount(uid, maxID); CHKERRQ(ierr);
  for (i = 0; i < 2*size; ++i) {
    ierr = UniqueIDGenerate(uid,&id); CHKERRQ(ierr);
    ierr = PetscSynchronizedPrintf(comm,"%d: %d\n", rank, id); CHKERRQ(ierr);
  }
  ierr = PetscSynchronizedFlush(comm); CHKERRQ(ierr);

  ierr = PetscPrintf(comm,"\nTesting Max ID\n"); CHKERRQ(ierr);
  if( rank == size-1 ) {
    ierr = UniqueIDSetStartCount(uid, PETSC_MAX_INT-4*size); CHKERRQ(ierr);
    for (i = 0; i < 5*size; ++i) {
//      ierr = UniqueIDGenerate(uid,&id); CHKERRQ(ierr);
      ierr = PetscSynchronizedPrintf(comm,"%d: %d  (%d)\n", rank, id, PETSC_MAX_INT-id); CHKERRQ(ierr);
    }
  }
  ierr = PetscSynchronizedFlush(comm); CHKERRQ(ierr);

  ierr = PetscPrintf(comm,"\nDestroying\n"); CHKERRQ(ierr);
  ierr = UniqueIDDestroy(uid); CHKERRQ(ierr);

	ierr = DCellFinalize(); CHKERRQ(ierr);
  return 0;
}
