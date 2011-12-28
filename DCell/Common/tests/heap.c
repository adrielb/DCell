#include "Common.h"

typedef struct {
  PetscReal phi;
  iCoor p;
} FMMNode;

static inline PetscReal MyComparator( void *parent, void *child );
void MyNodePrint( void *node );

int main(int argc, char **args)
{
  PetscErrorCode  ierr;
  ierr = DCellInit(); CHKERRQ(ierr);

  MemCache mc;
  MemCacheCreate(sizeof(FMMNode),1e5,&mc);

  Heap h;
  ierr = HeapCreate(MyComparator,&h); CHKERRQ(ierr);

  PetscRandom rnd;
  ierr = PetscRandomCreate(PETSC_COMM_SELF,&rnd); CHKERRQ(ierr);
  ierr = PetscRandomSetType(rnd,PETSCRAND48); CHKERRQ(ierr);

  int i;
  FMMNode *node;
  for ( i = 0; i < 10; ++i) {
    MemCacheAlloc(mc,&node);
    PetscRandomGetValue(rnd,&node->phi);
    ierr = HeapInsert(h,node); CHKERRQ(ierr);
    ierr = HeapCheck(h); CHKERRQ(ierr);
  }
  ierr = HeapPrint(h,MyNodePrint); CHKERRQ(ierr);
  for ( i = 0; i < 10; ++i) {
    MemCacheAlloc(mc,&node);
    PetscRandomGetValue(rnd,&node->phi);
    node->phi = -1*node->phi;
    ierr = HeapInsert(h,node); CHKERRQ(ierr);
    ierr = HeapCheck(h); CHKERRQ(ierr);
  }
  ierr = HeapPrint(h,MyNodePrint); CHKERRQ(ierr);

  node = 0;
  for ( i = 0; i < 20; ++i) {
    ierr = HeapPop(h,&node); CHKERRQ(ierr);
    printf("%d: %f\n",i,node->phi);
    ierr = MemCacheFree(mc,node); CHKERRQ(ierr);
    ierr = HeapCheck(h); CHKERRQ(ierr);
  }

  int t;
  int big = 100000;
  PetscLogDouble t1, tend;
  for ( t = 0; t < 10; ++t) {
    PetscGetTime(&t1);
    for ( i = 0; i < big; ++i) {
      MemCacheAlloc(mc,&node);
      PetscRandomGetValue(rnd,&node->phi);
      ierr = HeapInsert(h,node); CHKERRQ(ierr);
    }
    ierr = HeapCheck(h); CHKERRQ(ierr);
    for ( i = 0; i < big; ++i) {
      ierr = HeapPop(h,&node); CHKERRQ(ierr);
      ierr = MemCacheFree(mc,node); CHKERRQ(ierr);
    }
    PetscGetTime(&tend);
    printf("TIME: %f\n", tend-t1);
  }
  ierr = PetscRandomDestroy(rnd); CHKERRQ(ierr);
  ierr = HeapDestroy(h); CHKERRQ(ierr);
  ierr = MemCacheDestroy(mc); CHKERRQ(ierr);
  ierr = DCellFinalize(); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

inline PetscReal MyComparator( void *parent, void *child ) {
  return ((FMMNode*)parent)->phi - ((FMMNode*)child)->phi;
}

void MyNodePrint( void *node ) {
  FMMNode *n = (FMMNode*) node;
  printf(" %1.2f ",n->phi);
}
