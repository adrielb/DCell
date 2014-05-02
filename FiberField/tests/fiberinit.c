#include "FiberField.h"

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc, char **args) {

  PetscErrorCode ierr;
  ierr = DCellInit(); CHKERRQ(ierr);
  FiberField fibers;
  ierr = FiberFieldCreate( MPI_COMM_WORLD, &fibers); CHKERRQ(ierr);
  Vertex v0, v1, v2, v3;
  ierr = VertexCreate( fibers, &v0); CHKERRQ(ierr);
  ierr = VertexCreate( fibers, &v1); CHKERRQ(ierr);
  ierr = VertexCreate( fibers, &v2); CHKERRQ(ierr);
  ierr = VertexCreate( fibers, &v3); CHKERRQ(ierr);

  EdgeType e = 1;
  ierr = VertexAddEdge( v0, v1, e ); CHKERRQ(ierr);
  ierr = VertexAddEdge( v2, v3, e ); CHKERRQ(ierr);
 
  int i;
  int len = ArrayLength( fibers->fibers );
  printf("len: %d\n", len);
  struct _Vertex* vs = ArrayGetData(fibers->fibers);
  for (i = 0; i < len; i++) {
    printf("ID: %d\n", vs[i].id );
  }

  ierr = FiberFieldDestroy(fibers); CHKERRQ(ierr);
  ierr = DCellFinalize(); CHKERRQ(ierr);
  return 0;
}


