#include "FiberField.h"

#undef __FUNCT__
#define __FUNCT__ "FiberFieldCreate"
PetscErrorCode FiberFieldCreate(MPI_Comm comm, FiberField *fibers)
{
  FiberField f;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscNew(struct _FiberField, &f); CHKERRQ(ierr);
  ierr = ArrayCreate("verts",sizeof(struct _Vertex),&f->verts); CHKERRQ(ierr);
  ierr = ArrayCreate("edges",sizeof(struct _Vertex),&f->edges); CHKERRQ(ierr);
  ierr = UniqueIDCreate( &f->vid ); CHKERRQ(ierr);
  f->comm = comm;

  f->drag = 1;
  f->mass = 1;

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
  ierr = ArrayDestroy(fibers->verts); CHKERRQ(ierr);
  ierr = ArrayDestroy(fibers->edges); CHKERRQ(ierr);
  ierr = PetscFree( fibers ); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "FiberFieldPrint"
PetscErrorCode FiberFieldPrint( FiberField fibers )
{
  int i,j;
  int len = ArrayLength( fibers->verts );
  PetscErrorCode ierr;

  PetscFunctionBegin;
  struct _Vertex* vs = ArrayGetData(fibers->verts);
  for (i = 0; i < len; i++) {
    ierr = PetscSynchronizedPrintf( fibers->comm, "ID: %d; ", vs[i].vID ); CHKERRQ(ierr);
    for (j = 0; j < MAXEDGES; j++) {
      /*printf("%d ", vs[i].vID[j] ); */
      ierr = PetscSynchronizedPrintf( fibers->comm, "%d ", vs[i].eID[j] ); CHKERRQ(ierr);
    }
    ierr = PetscSynchronizedPrintf( fibers->comm, "\n" ); CHKERRQ(ierr);
  }
  ierr = PetscSynchronizedFlush(fibers->comm); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
