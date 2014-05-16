#include "FiberField.h"

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc, char **args) {
  PetscErrorCode ierr;
  ierr = DCellInit(); CHKERRQ(ierr);
  MPI_Comm comm = PETSC_COMM_WORLD;
  FiberField fibers;
  ierr = FiberFieldCreate( comm, &fibers); CHKERRQ(ierr);


  MPI_Datatype verttype; 

#define LEN 4
  struct _VertexMPI v;
  MPI_Aint disp[LEN];
  MPI_Get_address( &v     , disp);   // vID
  MPI_Get_address( &v.eID , disp+1); // eID
  MPI_Get_address( &v.x   , disp+2); // Coor x
  MPI_Get_address( &v.v   , disp+3); // Coor v

  MPI_Datatype types[] = {
    MPI_INT,    // vID
    MPI_INT,    // eID
    MPI_DOUBLE, // Coor x
    MPI_DOUBLE  // Coor v
  };
  int blocklen[] = { 
    1,        // vID
    MAXEDGES, // eID
    3,        // Coor x
    3         // Coor v
  };
  
  MPI_Type_create_struct(LEN, blocklen, disp, types, &verttype );

  MPI_Type_commit( &verttype );

  
  MPI_Type_free( &verttype ); 
  ierr = FiberFieldDestroy(fibers); CHKERRQ(ierr);
  ierr = DCellFinalize(); CHKERRQ(ierr);
  return 0;
}
