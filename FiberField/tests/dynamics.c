#include "FiberField.h"

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc, char **args) {

  PetscErrorCode ierr;
  ierr = DCellInit(); CHKERRQ(ierr);
  FiberField fibers;
  ierr = FiberFieldCreate( PETSC_COMM_WORLD, &fibers); CHKERRQ(ierr);

  //1. add fibers to system
  //2. set them in motion

  Coor gmin = {-1,-1,-1};
  Coor gmax = {1,1,1};
  fibers->dh = 0.1;
  fibers->globalBounds = (BoundingBox){.min = gmin, .max = gmax };
  ierr = FiberFieldSetup(fibers); CHKERRQ(ierr);

  int rank;
  MPI_Comm_rank( PETSC_COMM_WORLD, &rank );
  
  if( rank == 0 ) {

    Vertex v0;
    ierr = FiberFieldAddVertex( fibers, &v0); CHKERRQ(ierr);
    v0->X.x = 0;
    v0->X.y = 0;
    v0->X.z = 0;
    v0->V.x = 0;
    v0->V.y = 0;
    v0->V.z = 0;

    Vertex v1;
    ierr = FiberFieldAddVertex( fibers, &v1); CHKERRQ(ierr);
    v1->X.x = 0;
    v1->X.y = 1;
    v1->X.z = 0;
    v1->V.x = 0;
    v1->V.y = 0;
    v1->V.z = 0;

    Vertex v2;
    ierr = FiberFieldAddVertex( fibers, &v2); CHKERRQ(ierr);
    v2->X.x = 1;
    v2->X.y = 1;
    v2->X.z = 0;
    v2->V.x = 0;
    v2->V.y = 0;
    v2->V.z = 0;
  
    EdgeType collagen = 1;
    PetscReal l0 = 3.1;

    ierr = FiberFieldAddEdge( fibers, v0, v1, collagen, l0 ); CHKERRQ(ierr);
    ierr = FiberFieldAddEdge( fibers, v1, v2, collagen, l0 ); CHKERRQ(ierr);

  } 


  ierr = FiberFieldPrint( fibers ); CHKERRQ(ierr);

  fibers->dt = .1;

  int t;
  int timax = 10;
  for (t = 0; t < timax; t++) {
    ierr = FiberFieldSolve(fibers); CHKERRQ(ierr);
    /*ierr = FiberFieldView( fibers ); CHKERRQ(ierr);*/
    ierr = FiberFieldWrite( fibers, t); CHKERRQ(ierr);
  }
 
  ierr = FiberFieldDestroy(fibers); CHKERRQ(ierr);
  ierr = DCellFinalize(); CHKERRQ(ierr);
  return 0;
}
