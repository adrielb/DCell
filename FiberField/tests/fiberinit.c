#include "FiberField.h"

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc, char **args) {

  PetscErrorCode ierr;
  ierr = DCellInit(); CHKERRQ(ierr);
  FiberField fibers;
  ierr = FiberFieldCreate( MPI_COMM_WORLD, &fibers); CHKERRQ(ierr);

  Coor gmin = {-1,-1,-1};
  Coor gmax = {1,1,1};
  fibers->dh = 0.1;
  fibers->globalBounds = (BoundingBox){.min = gmin, .max = gmax };

  FiberTypeID collagenVert;
  FiberTypeID collagenEdge;
  FiberTypeID collagenBendEdge;
  ierr = FiberFieldAddType( fibers,     "collagen vert",PETSC_FALSE,&collagenVert); CHKERRQ(ierr);
  ierr = FiberFieldAddType( fibers,     "collagen edge", PETSC_TRUE,&collagenEdge); CHKERRQ(ierr);
  ierr = FiberFieldAddType( fibers,"collagen bend edge", PETSC_TRUE,&collagenBendEdge); CHKERRQ(ierr);

  ierr = FiberFieldSetup(fibers); CHKERRQ(ierr);

  PetscReal l0 = 0.1;
  const int MAX_VERTS = 100;
  ierr = FiberFieldInitLocalFiber(fibers, MAX_VERTS, l0,
      collagenVert, collagenEdge, collagenBendEdge); CHKERRQ(ierr);
  ierr = FiberFieldWrite(fibers, 0); CHKERRQ(ierr);

  ierr = FiberFieldDestroy(fibers); CHKERRQ(ierr);
  ierr = DCellFinalize(); CHKERRQ(ierr);
  return 0;
}


