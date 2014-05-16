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

  //1. add fibers to system
  //2. set them in motion

  printf("Adding edge\n");
  EdgeType e = 1;
  ierr = VertexAddEdge( v0, v1, e ); CHKERRQ(ierr);
  ierr = VertexAddEdge( v2, v3, e ); CHKERRQ(ierr);
  ierr = FiberFieldPrint( fibers ); CHKERRQ(ierr);

  ierr = FiberField_Setup(fibers); CHKERRQ(ierr);
  ierr = FiberField_Step(fibers); CHKERRQ(ierr);
 

  ierr = FiberFieldDestroy(fibers); CHKERRQ(ierr);
  ierr = DCellFinalize(); CHKERRQ(ierr);
  return 0;
}
