#include "FiberField.h"

#undef __FUNCT__
#define __FUNCT__ "FiberFieldRedistribute"
PetscErrorCode FiberFieldRedistribute(FiberField fiberfield)
{
  const int numNei = is2D ? 9 : 27; // or 0 if uniproc
  int myrank, rank;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  MPI_Comm_rank(world->comm,&myrank);
  for each Vertex v
  {
    rank = VertexRank( v );
    if( rank != myrank ) {
      for( r in 27 )
      {
        if( nei[r] == rank )
        {
          AddVertex( idx[r], v);
        }
      }
    }
  }

  for ( r in 27 ) {
    MPI_Type_indexed(count,blocklens[],offsets[],old_type,&new_type);
    MPI_Type_commit(new_type);
    MPI_Send(vertices,1,new_type,nei[r],)

    MPI_Probe( source, tag, comm, &status);
    MPI_Get_Count( status, datatype, &count );
    ierr = ArraySetSize( particleBuffer, count ); CHKERRQ(ierr);
    MPI_Recv(particleBuffer,MAXELEM,particleType,source,tag,PETSC_COMM_WORLD);
    MPI_Type_free(new_type);
  }

  // OR: post all 27 sends simultaneously (non-blocking), then post 27 recv
  PetscFunctionReturn(0);
}

