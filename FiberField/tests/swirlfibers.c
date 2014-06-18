#include "FiberField.h"

void VelocityField( Coor X, Coor *V);

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc, char **args) {
  PetscErrorCode ierr;
  ierr = DCellInit(); CHKERRQ(ierr);
  FiberField fibers;
  MPI_Comm comm = PETSC_COMM_WORLD;
  ierr = FiberFieldCreate( comm, &fibers); CHKERRQ(ierr);


  // create world dims
  // add fibers
  // advect 
  // rebalance
  
  Coor gmin = {-1,-1,-1};
  Coor gmax = {1,1,1};

  fibers->dh = 0.1;
  fibers->globalBounds = (BoundingBox){.min = gmin, .max = gmax };
  ierr = FiberFieldSetup(fibers); CHKERRQ(ierr);

  BoundingBox bbox = fibers->localBounds;

  FiberTypeID collagenVert;
  FiberTypeID collagenEdge;
  FiberTypeID collagenBendEdge;
  ierr = FiberFieldAddType( fibers,     "collagen vert",PETSC_FALSE,&collagenVert) ; CHKERRQ(ierr);
  ierr = FiberFieldAddType( fibers,     "collagen edge", PETSC_TRUE,&collagenEdge) ; CHKERRQ(ierr);
  ierr = FiberFieldAddType( fibers,"collagen bend edge", PETSC_TRUE,&collagenBendEdge) ; CHKERRQ(ierr);

//init verts
  Array fiber;
  ierr = ArrayCreate("fiber", sizeof(Vertex), &fiber); CHKERRQ(ierr);
  PetscReal a; 
  Vertex v1;
  Coor d = { 
    bbox.max.x - bbox.min.x, 
    bbox.max.y - bbox.min.y, 
    bbox.max.z - bbox.min.z 
  };
  const int flen = 80;
  d.x = d.x / (flen+1);
  d.y = d.y / (flen+1);
  d.z = d.z / (flen+1);

  ierr = PetscInfo3(0, "d = {%f,%f,%f}\n", d.x, d.y, d.z); CHKERRQ(ierr);

  for (a = 0; a < flen; a++ ) {
    ierr = FiberFieldAddVertex( fibers, collagenVert, &v1); CHKERRQ(ierr);
    v1->X.x = d.x * a + bbox.min.x + d.x;
    v1->X.y = d.y * a + bbox.min.y + d.y;
    v1->X.z = d.z * a + bbox.min.z + d.z;

    ierr = ArrayAppendPtr( fiber, v1); CHKERRQ(ierr);
  }

  const PetscReal norm_d = PetscSqrtReal(d.x*d.x + d.y*d.y + d.z*d.z);
  const PetscReal l0 = 0.9*norm_d;
  int i;
  const int len = ArrayLength( fiber );
  Vertex *v = ArrayGetData( fiber );
  for (i = 0; i < len-1; i++) {
    ierr = FiberFieldAddEdge( fibers, v[i], v[i+1], collagenEdge, l0 ); CHKERRQ(ierr);
  }
  for (i = 1; i < len-1; i++) {
    ierr = FiberFieldAddEdge( fibers, v[i-1], v[i+1], collagenBendEdge, 2*l0 ); CHKERRQ(ierr);
  }

// advect verts
  int ti;
  int timax = 1500;
  fibers->dt = 0.01;
  fibers->fluidDrag = 10;
  ierr = FiberFieldSetFluidVelocityEvaluator( fibers, FiberField_CircularFluidVelocity); CHKERRQ(ierr);

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
