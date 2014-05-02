#include "FiberField.h"

#undef __FUNCT__
#define __FUNCT__ "FiberFieldCreate"
PetscErrorCode FiberFieldCreate(MPI_Comm comm, FiberField *fibers)
{
  FiberField f;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscNew(struct _FiberField, &f); CHKERRQ(ierr);
  ierr = ArrayCreate("fibers",sizeof(struct _Vertex),&f->fibers); CHKERRQ(ierr);
  ierr = UniqueIDCreate( &f->vid ); CHKERRQ(ierr);

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
  ierr = ArrayDestroy(fibers->fibers); CHKERRQ(ierr);
  ierr = PetscFree( fibers ); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "FiberFieldPrint"
PetscErrorCode FiberFieldPrint( FiberField fibers )
{
  int i,j;
  int len = ArrayLength( fibers->fibers );
  //PetscErrorCode ierr;

  PetscFunctionBegin;
  struct _Vertex* vs = ArrayGetData(fibers->fibers);
  for (i = 0; i < len; i++) {
    printf("ID: %d; ", vs[i].ID );
    for (j = 0; j < MAXEDGES; j++) {
      printf("%d ", vs[i].vID[j] ); 
    }
    printf( "\n" );
  }
  PetscFunctionReturn(0);
}
