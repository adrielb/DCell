#include "LevelSetMethod.h"

// Common workspace shared by all level set objects

// Min-heap
static inline PetscReal FMMComparator( void *parent, void *child ) {
  return ((FMMNode)parent)->phi - ((FMMNode)child)->phi;
}

typedef struct _LevelSetCommon {
  MemCache mc;
  Heap heap;
  Grid tmp;
  Grid velgrid;
  Array val, coor, idx;
} *LevelSetCommon;

const

PetscErrorCode LevelSetInitialize( )
{
  ierr = MemCacheCreate(sizeof(struct _FMMNode),1024,&mc); CHKERRQ(ierr);
  ierr = HeapCreate(FMMComparator, &heap); CHKERRQ(ierr);
  ierr = GridCreate(ls->tmp->n,1,&tmp); CHKERRQ(ierr);
  ierr = GridCreate(ls->phi->n,dof,&velgrid); CHKERRQ(ierr);
  ierr = ArrayCreate(sizeof(PetscReal),len,&val); CHKERRQ(ierr);
  ierr = ArrayCreate(sizeof(int*)     ,len,&idx); CHKERRQ(ierr);
  ierr = ArrayCreate(sizeof(int)*4    ,len,&coor); CHKERRQ(ierr); // [z,y,x,d]
}

PetscErrorCode LevelSetFinalize()
{
  PetscRegisterFinalize(LevelSetFinalize)
}
