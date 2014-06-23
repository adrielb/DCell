#include "FiberField.h"

int main(int argc, char **args)
{
  FiberField fibers;
  PetscErrorCode ierr;
  ierr = DCellInit(); CHKERRQ(ierr);
  ierr = FiberFieldCreate( MPI_COMM_WORLD, &fibers); CHKERRQ(ierr);

  Coor gmin = {-1,-1,-1};
  Coor gmax = {1,1,1};
  fibers->dh = 0.1;
  fibers->globalBounds = (BoundingBox){.min = gmin, .max = gmax };

  FiberTypeID collagenVert;
  FiberTypeID collagenEdge;
  FiberTypeID collagenBendEdge;
  ierr = FiberFieldAddType( fibers,     "collagen vert",PETSC_FALSE,&collagenVert) ; CHKERRQ(ierr);
  ierr = FiberFieldAddType( fibers,     "collagen edge", PETSC_TRUE,&collagenEdge) ; CHKERRQ(ierr);
  ierr = FiberFieldAddType( fibers,"collagen bend edge", PETSC_TRUE,&collagenBendEdge) ; CHKERRQ(ierr);

  ierr = FiberFieldSetup(fibers); CHKERRQ(ierr);

  int i;
  PetscReal l0 = 0.1;
  int MAX_VERTS = 100;
  int NUM_FIBERS_PER_PROC = 50;
  ierr = PetscOptionsGetInt(0,"-num_fibers_per_proc",&NUM_FIBERS_PER_PROC ,0); CHKERRQ(ierr);
  ierr = PetscOptionsGetInt(0,"-max_verts",&MAX_VERTS,0); CHKERRQ(ierr);
  for (i = 0; i < NUM_FIBERS_PER_PROC; i++) {
    ierr = FiberFieldInitLocalFiber(fibers, MAX_VERTS, l0,
        collagenVert, collagenEdge, collagenBendEdge); CHKERRQ(ierr);
  }

// advect verts
  int ti;
  int timax = 800;
  fibers->dt = 0.5;
  fibers->fluidDrag = 10;
  /*ierr = FiberFieldSetFluidVelocityEvaluator( fibers, FiberField_CircularFluidVelocity); CHKERRQ(ierr);*/

  ierr = PetscOptionsGetInt(0,"-timax",&timax ,0); CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(0,"-fiber_dt",&fibers->dt,0); CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(0,"-fluidDrag",&fibers->fluidDrag,0); CHKERRQ(ierr);

  ierr = FiberFieldWrite( fibers, 0 ); CHKERRQ(ierr);

  for (ti = 1; ti < timax; ti++) {
    ierr = PetscInfo1(0, "time = %d\n", ti); CHKERRQ(ierr);
    ierr = FiberFieldSolve(fibers); CHKERRQ(ierr);
    ierr = FiberFieldWrite( fibers, ti ); CHKERRQ(ierr);
  }

  ierr = FiberFieldDestroy(fibers); CHKERRQ(ierr);
  ierr = DCellFinalize(); CHKERRQ(ierr);
  return 0;
}
