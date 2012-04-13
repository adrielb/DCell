#include "Common.h"

struct _UniqueID {
  UniqueIDType localCount;
  UniqueIDType MAX_ID;
  int rank;
  int size;
};

#undef __FUNCT__
#define __FUNCT__ "UniqueIDCreate"
PetscErrorCode UniqueIDCreate( UniqueID *uid )
{
  UniqueID uniqueID;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscNew(struct _UniqueID, &uniqueID); CHKERRQ(ierr);
  ierr = MPI_Comm_size(PETSC_COMM_WORLD,&uniqueID->size); CHKERRQ(ierr);
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&uniqueID->rank); CHKERRQ(ierr);
  uniqueID->localCount = 0;
  uniqueID->MAX_ID = PETSC_MAX_INT - 2 * uniqueID->size;

  *uid = uniqueID;
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "UniqueIDDestroy"
PetscErrorCode UniqueIDDestroy( UniqueID uid )
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscFree(uid); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "UniqueIDSetStartCount"
PetscErrorCode UniqueIDSetStartCount( UniqueID uid, UniqueIDType maxID )
{
  PetscFunctionBegin;
  uid->localCount = (maxID - uid->rank) / uid->size + 1;
  PetscSynchronizedPrintf(MPI_COMM_WORLD,"lc: %d\n", uid->localCount);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "UniqueIDGenerate"
PetscErrorCode UniqueIDGenerate( UniqueID uid, UniqueIDType *id )
{
  PetscFunctionBegin;

  *id = uid->rank + uid->size * uid->localCount;
  uid->localCount++;

  if( *id >= uid->MAX_ID ) {
    SETERRABORT(PETSC_COMM_WORLD,0,"UniqueID: Max ID reached");
  }

  PetscFunctionReturn(0);
}

