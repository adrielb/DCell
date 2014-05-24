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
  // add verticies
  // advect 
  // rebalance
  
  Coor gmin = {-1,-1,-1};
  Coor gmax = {1,1,1};

  fibers->dh = 0.1;
  fibers->globalBounds = (BoundingBox){.min = gmin, .max = gmax };
  ierr = FiberFieldSetup(fibers); CHKERRQ(ierr);
  gmin.x += PETSC_SMALL;
  gmin.y += PETSC_SMALL;
  gmin.z += PETSC_SMALL;
  gmax.x -= PETSC_SMALL;
  gmax.y -= PETSC_SMALL;
  gmax.z -= PETSC_SMALL;

  BoundingBox bbox = fibers->localBounds;

//init verts
  PetscReal a,b;
  Vertex v;
  for (a  = 0; a < 1; a+=0.1) {
    for (b = 0; b < 1; b+=0.1) {
      ierr = VertexCreate( fibers, &v); CHKERRQ(ierr);
      v->X.x = a*bbox.max.x + (1-a)*bbox.min.x;
      v->X.y = (bbox.max.y + bbox.min.y) / 2;
      v->X.z = b*bbox.max.z + (1-b)*bbox.min.z;
    }
  }

// advect verts
  int i;
  int ti;
  int timax = 300;
  PetscReal dt = 0.05;
  for (ti = 0; ti < timax; ti++) {
    ierr = PetscInfo1(0, "time = %d\n", ti); CHKERRQ(ierr);

    for (i = 0; i < ArrayLength(fibers->verts); i++) {
      ierr = ArrayGet( fibers->verts, i, &v ); CHKERRQ(ierr);

      VelocityField( v->X, &v->V );

      v->X.x = v->X.x + dt * v->V.x;
      v->X.y = v->X.y + dt * v->V.y;
      v->X.z = v->X.z + dt * v->V.z;

      v->X.x = v->X.x > gmax.x ? gmax.x : v->X.x;
      v->X.y = v->X.y > gmax.y ? gmax.y : v->X.y;
      v->X.z = v->X.z > gmax.z ? gmax.z : v->X.z;
      v->X.x = v->X.x < gmin.x ? gmin.x : v->X.x;
      v->X.y = v->X.y < gmin.y ? gmin.y : v->X.y;
      v->X.z = v->X.z < gmin.z ? gmin.z : v->X.z;

      ierr = PetscInfo1(0, "vID = %d\n", v->vID ); CHKERRQ(ierr);
      ierr = PetscInfo3(0, "X = %f,%f,%f\n", v->X.x,v->X.y,v->X.z ); CHKERRQ(ierr);
    }


    ierr = FiberField_SpatiallyBalance( fibers ); CHKERRQ(ierr);

  }

  ierr = FiberFieldDestroy(fibers); CHKERRQ(ierr);
  ierr = DCellFinalize(); CHKERRQ(ierr);
  return 0;
}

void VelocityField( Coor X, Coor *V) {
  V->x = -X.z;
  V->y =  0;
  V->z =  X.x; 
}
