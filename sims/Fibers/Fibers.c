#include "FiberField.h"

int main(int argc, char **args)
{
  FiberField fibers;
  PetscErrorCode ierr;
  ierr = DCellInit(); CHKERRQ(ierr);
  ierr = FiberFieldCreate( MPI_COMM_WORLD, &fibers); CHKERRQ(ierr);

  Vertex v0, v1, v2, v3;
  ierr = VertexCreate( fibers, &v0); CHKERRQ(ierr);
  ierr = VertexCreate( fibers, &v1); CHKERRQ(ierr);
  ierr = VertexCreate( fibers, &v2); CHKERRQ(ierr);
  ierr = VertexCreate( fibers, &v3); CHKERRQ(ierr);

  EdgeType e = 1;
  ierr = VertexAddEdge( v0, v1, e ); CHKERRQ(ierr);
  ierr = VertexAddEdge( v2, v3, e ); CHKERRQ(ierr);

  ierr = FiberFieldDestroy(fibers); CHKERRQ(ierr);
  ierr = DCellFinalize(); CHKERRQ(ierr);
  return 0;
}
