#include "FiberField.h"

#undef __FUNCT__
#define __FUNCT__ "FiberFieldCreate"
PetscErrorCode FiberFieldCreate(MPI_Comm comm, FiberField *fibers)
{
  FiberField f;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscNew(struct _FiberField, &f); CHKERRQ(ierr);
  ierr = ArrayCreate("fibers",sizeof(Vertex),1000,&f->fibers); CHKERRQ(ierr);
  ierr = MemCacheCreate(sizeof(struct _Edge),1000,&f->mcEdges); CHKERRQ(ierr);
  ierr = MemCacheCreate(sizeof(struct _Vertex),1000,&f->mcVerticies); CHKERRQ(ierr);
  ierr = UniqueIDCreate( &f->vid ); CHKERRQ(ierr);
  ierr = UniqueIDCreate( &f->eid ); CHKERRQ(ierr);

  *fibers = f;
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "FiberFieldDestroy"
PetscErrorCode FiberFieldDestroy(FiberField fibers)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = UniqueIDDestroy(fibers->vid); CHKERRQ(ierr);
  ierr = UniqueIDDestroy(fibers->eid); CHKERRQ(ierr);
  ierr = MemCacheDestroy(fibers->mcEdges); CHKERRQ(ierr);
  ierr = MemCacheDestroy(fibers->mcVerticies); CHKERRQ(ierr);
  ierr = ArrayDestroy(fibers->fibers); CHKERRQ(ierr);
  ierr = PetscFree( fibers ); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
