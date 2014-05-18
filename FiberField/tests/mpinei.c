#include "FiberField.h"

//#define PRINT_NEI_RANKS
#define PRINT_BBOX

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc, char **args) {
  PetscErrorCode ierr;
  ierr = DCellInit(); CHKERRQ(ierr);

  int i;
  UNUSED(i);
  FiberField fibers;
  MPI_Comm comm = PETSC_COMM_WORLD;
  ierr = FiberFieldCreate( comm, &fibers); CHKERRQ(ierr);

  int size;
  MPI_Comm_size( comm, &size);
  ierr = PetscPrintf( comm, "Comm size: %d\n", size); CHKERRQ(ierr);
  int rank; 
  MPI_Comm_rank( comm, &rank ); 

  PetscReal dx = 0.1;
  Coor dX = {dx, dx, dx};
  Coor len = {5,3,4};
  iCoor dims = { 
    len.x/dX.x,
    len.y/dX.y, 
    len.z/dX.z
  };
  const int dof = 1;
  DM da;

  ierr = DMDACreate3d(comm,
      DMDA_BOUNDARY_NONE,DMDA_BOUNDARY_NONE,DMDA_BOUNDARY_NONE,
      DMDA_STENCIL_BOX,
      dims.x,dims.y,dims.z, 
      PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE,
      dof,1,  0,0,0, &da); CHKERRQ(ierr);

  const PetscMPIInt *nei;
  ierr = DMDAGetNeighbors(da, &nei); CHKERRQ(ierr);

#ifdef PRINT_NEI_RANKS
  ierr = PetscSynchronizedPrintf( comm, "[%d] ",rank); CHKERRQ(ierr);
  for (i = 0; i < 27; i++) {
    if (nei[i] == -2 ) {
      ierr = PetscSynchronizedPrintf( comm, "  "); CHKERRQ(ierr);
    } else { 
      ierr = PetscSynchronizedPrintf( comm, "%d ", nei[i]); CHKERRQ(ierr);
    }    
  }
  ierr = PetscSynchronizedPrintf( comm, "\n"); CHKERRQ(ierr);
  ierr = PetscSynchronizedFlush( comm); CHKERRQ(ierr);
#endif

  DMDALocalInfo info;
  DMDAGetLocalInfo(da, &info);
  Coor lmin = {
    dX.x * info.xs,
    dX.y * info.ys,
    dX.z * info.zs
  };
  Coor lmax = {
    dX.x * (info.xs + info.xm),
    dX.y * (info.ys + info.ym),
    dX.z * (info.zs + info.zm)
  };
  
#ifdef PRINT_BBOX
  ierr = PetscSynchronizedPrintf( comm, "[%d] %d, %d, %d\n", rank, info.xs,info.ys,info.zs ); CHKERRQ(ierr);
  ierr = PetscSynchronizedFlush( comm); CHKERRQ(ierr);
  ierr = PetscSynchronizedPrintf( comm, "[%d] %d, %d, %d\n", rank, info.xm,info.ym,info.zm ); CHKERRQ(ierr);
  ierr = PetscSynchronizedFlush( comm); CHKERRQ(ierr);
  ierr = PetscSynchronizedPrintf( comm, "[%d] %f, %f, %f\n", rank, lmin.x,lmin.y,lmin.z ); CHKERRQ(ierr);
  ierr = PetscSynchronizedPrintf( comm, "[%d] %f, %f, %f\n", rank, lmax.x,lmax.y,lmax.z ); CHKERRQ(ierr);
  ierr = PetscSynchronizedFlush( comm); CHKERRQ(ierr);
#endif

  UNUSED(lmax);
  Coor p = {3.7,1,3};
  iCoor x = {
    floor((p.x - lmin.x) / (lmax.x-lmin.x))+1, 
    floor((p.y - lmin.y) / (lmax.y-lmin.y))+1, 
    floor((p.z - lmin.z) / (lmax.z-lmin.z))+1 
  };
  ierr = PetscSynchronizedPrintf( comm, "[%d] %d, %d, %d  ", rank, x.x, x.y, x.z ); CHKERRQ(ierr);
  if (x.x < 0 || x.x > 2 ||
      x.y < 0 || x.y > 2 ||
      x.z < 0 || x.z > 2 ) {
    ierr = PetscSynchronizedPrintf( comm, "X\n" ); CHKERRQ(ierr);
  } else { 
    ierr = PetscSynchronizedPrintf( comm, "%d\n", nei[x.x + 3*x.y + 9*x.z] ); CHKERRQ(ierr);
  }
  ierr = PetscSynchronizedFlush( comm); CHKERRQ(ierr);
  


  //PetscReal lmin[3], lmax[3];
  //ierr = DMDAGetBoundingBox(da, gmin, gmax); CHKERRQ(ierr);
  //ierr = PetscPrintf(comm, "gmin: %d, %d, %d\n" , gmin[0],gmin[1],gmin[2]); CHKERRQ(ierr);
  //ierr = PetscPrintf(comm, "gmax: %d, %d, %d\n" , gmax[0],gmax[1],gmax[2]); CHKERRQ(ierr);

//DMDAGetLocalBoundingBox

  ierr = FiberFieldDestroy(fibers); CHKERRQ(ierr);
  ierr = DCellFinalize(); CHKERRQ(ierr);
  return 0;
}
