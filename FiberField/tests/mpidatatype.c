#include "FiberField.h"

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc, char **args) {
  PetscErrorCode ierr;
  ierr = DCellInit(); CHKERRQ(ierr);
  MPI_Comm comm = PETSC_COMM_WORLD;
  FiberField fibers;
  ierr = FiberFieldCreate( comm, &fibers); CHKERRQ(ierr);

  int i;
  MPI_Datatype verttype; 

#define LEN 4
  struct _VertexMPI v;
  MPI_Aint disp[LEN];
  MPI_Get_address( &v     , disp);   // vID
  MPI_Get_address( &v.eID , disp+1); // eID
  MPI_Get_address( &v.X   , disp+2); // Coor x
  MPI_Get_address( &v.V   , disp+3); // Coor v

  disp[0] = offsetof(struct _VertexMPI, vID);
  disp[1] = offsetof(struct _VertexMPI, eID);
  disp[2] = offsetof(struct _VertexMPI, X);
  disp[3] = offsetof(struct _VertexMPI, V);

  for (i = 0; i < LEN; i++) {
    /*disp[i] = disp[i] - disp[0];*/
    printf("disp[%d] = %d\n", i, (int)disp[i] );
  }

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

  int tag = 0;
  int rank;
  MPI_Comm_rank( PETSC_COMM_WORLD, &rank );
  Array buf;
  ArrayCreate( "buf", sizeof( struct _VertexMPI ), &buf );
  ArraySetSize( buf, 1);
  if (rank == 0) {
    ierr = MPI_Send(ArrayGetData(buf), ArrayLength(buf), verttype, 1, tag, PETSC_COMM_WORLD); CHKERRQ(ierr);
  }
  if (rank == 1) {
    ierr = MPI_Recv(ArrayGetData(buf), ArrayLength(buf), verttype, 0, tag, PETSC_COMM_WORLD, MPI_STATUS_IGNORE ); CHKERRQ(ierr);
    printf("\n");
  }
  

  MPI_Type_free( &verttype ); 
  ierr = FiberFieldDestroy(fibers); CHKERRQ(ierr);
  ierr = DCellFinalize(); CHKERRQ(ierr);
  return 0;
}
