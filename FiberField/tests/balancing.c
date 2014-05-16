#include "FiberField.h"

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc, char **args) {
  PetscErrorCode ierr;
  ierr = DCellInit(); CHKERRQ(ierr);
  MPI_Comm comm = PETSC_COMM_WORLD;
  FiberField fibers;
  ierr = FiberFieldCreate( comm, &fibers); CHKERRQ(ierr);

  // Randomly seed verts 
  // put into motion
  // rebalance 

  int i;
  int nLocal = 3;
  Vertex v;
  for (i = 0; i < nLocal; i++) {
    ierr = VertexCreate( fibers, &v ); CHKERRQ(ierr);
  }

  /*ierr = FiberFieldPrint( fibers ); CHKERRQ(ierr);*/

#define LEN 6
  int rank;
  int size;
  MPI_Comm_rank( comm, &rank);
  MPI_Comm_size( comm, &size);
  int sendbuf[LEN];
  int recvbuf[LEN];

  for (i = 0; i < LEN; i++) {
    sendbuf[i] = i*size+rank;
    sendbuf[i] = rank;
    recvbuf[i] = -1;
  }
  
  ierr = PetscPrintf( comm, "comm size: %d\n" , size); CHKERRQ(ierr);
  
/********************************************************************************/
  ierr = PetscPrintf( comm, "MPI All to all\n" , size); CHKERRQ(ierr);
  ierr = PetscSynchronizedPrintf( comm, "[%d] sendbuf: ", rank); CHKERRQ(ierr);
  for (i = 0; i < LEN; i++) {
    ierr = PetscSynchronizedPrintf( comm, "%d ", sendbuf[i]); CHKERRQ(ierr);
  }
  ierr = PetscSynchronizedPrintf( comm, "\n"); CHKERRQ(ierr);
  ierr = PetscSynchronizedFlush(comm); CHKERRQ(ierr);

  ierr = MPI_Alltoall(sendbuf, 2, MPI_INT,
                      recvbuf, 2, MPI_INT, 
                      PETSC_COMM_WORLD); CHKERRQ(ierr);

  ierr = PetscSynchronizedPrintf( comm, "[%d] recvbuf: ", rank); CHKERRQ(ierr);
  for (i = 0; i < LEN; i++) {
    ierr = PetscSynchronizedPrintf( comm, "%d ", recvbuf[i]); CHKERRQ(ierr);
  }
  ierr = PetscSynchronizedPrintf( comm, "\n"); CHKERRQ(ierr);
  ierr = PetscSynchronizedFlush(comm); CHKERRQ(ierr);

/********************************************************************************/
  ierr = PetscPrintf( comm, "MPI All to all variable size\n" , size); CHKERRQ(ierr);
  ierr = PetscSynchronizedPrintf( comm, "[%d] sendbuf: ", rank); CHKERRQ(ierr);
  for (i = 0; i < LEN; i++) {
    ierr = PetscSynchronizedPrintf( comm, "%d ", sendbuf[i]); CHKERRQ(ierr);
  }
  ierr = PetscSynchronizedPrintf( comm, "\n"); CHKERRQ(ierr);
  ierr = PetscSynchronizedFlush(comm); CHKERRQ(ierr);

  int sendcounts[LEN];
  int recvcounts[LEN];
  int sdispls[3] = {0,2,4};
  int rdispls[3] = {0,2,4};
  for (i = 0; i < LEN; i++) {
    sendcounts[i] = 2;
    recvcounts[i] = 2;
  }
  
  ierr = MPI_Alltoallv(sendbuf, sendcounts, sdispls, MPI_INT,
                       recvbuf, recvcounts, rdispls, MPI_INT, 
                       PETSC_COMM_WORLD); CHKERRQ(ierr);

  ierr = PetscSynchronizedPrintf( comm, "[%d] recvbuf: ", rank); CHKERRQ(ierr);
  for (i = 0; i < LEN; i++) {
    ierr = PetscSynchronizedPrintf( comm, "%d ", recvbuf[i]); CHKERRQ(ierr);
  }
  ierr = PetscSynchronizedPrintf( comm, "\n"); CHKERRQ(ierr);
  ierr = PetscSynchronizedFlush(comm); CHKERRQ(ierr);
  
/********************************************************************************/
  ierr = PetscPrintf( comm, "MPI point to point\n" , size); CHKERRQ(ierr);

  int tag = 0;
  int nei[2];
  MPI_Request reqSend[2];
  MPI_Status status;

  if (rank == 0) { nei[0] = 1; nei[1] = 2; }
  if (rank == 1) { nei[0] = 0; nei[1] = 2; }
  if (rank == 2) { nei[0] = 0; nei[1] = 1; }

  Array sendbufs[2];
  ierr = ArrayCreate("send0", sizeof(int), &sendbufs[0]); CHKERRQ(ierr);
  ierr = ArrayCreate("send1", sizeof(int), &sendbufs[1]); CHKERRQ(ierr);
  for (i = 0; i < 19; i++) {
    int *elem;
    ierr = ArrayAppend(sendbufs[0], &elem); CHKERRQ(ierr);
    *elem = rank;
  }
  for (i = 0; i < 9; i++) {
    int *elem;
    ierr = ArrayAppend(sendbufs[1], &elem); CHKERRQ(ierr);
    *elem = rank;
  }

  int numNei = 2;
  for (i = 0; i < numNei; i++) {
    // send 
    ierr = MPI_Isend(ArrayGetData(sendbufs[i]), ArrayLength(sendbufs[i]), MPI_INT, nei[i], tag, comm, &reqSend[i] ); CHKERRQ(ierr);
    /*ierr = PetscSynchronizedPrintf( comm, "[%d] Sending count: %d\n", rank, ArrayLength(sendbufs[i]) ); CHKERRQ(ierr);*/
  }
  /*ierr = PetscSynchronizedFlush(comm); CHKERRQ(ierr);*/


  MPI_Request reqRecv[2];
  Array recvbufs[2];
  ierr = ArrayCreate("recv0", sizeof(int), &recvbufs[0]); CHKERRQ(ierr);
  ierr = ArrayCreate("recv1", sizeof(int), &recvbufs[1]); CHKERRQ(ierr);
  int count;
  for (i = 0; i < numNei; i++) {
    ierr = MPI_Probe(MPI_ANY_SOURCE, tag, comm, &status); CHKERRQ(ierr);
    ierr = MPI_Get_count( &status, MPI_INT, &count); CHKERRQ(ierr);
    ierr = ArraySetSize( recvbufs[i], count); CHKERRQ(ierr);
    /*ierr = PetscSynchronizedPrintf( comm, "[%d] Recv count: %d\n", rank, ArrayLength(recvbufs[i]) ); CHKERRQ(ierr);*/
    ierr = MPI_Irecv( ArrayGetData(recvbufs[i]), count, MPI_INT, status.MPI_SOURCE, tag, comm, &reqRecv[i]); CHKERRQ(ierr);
  }
  /*ierr = PetscSynchronizedFlush(comm); CHKERRQ(ierr);*/

  MPI_Status statusRecv[2];
  MPI_Status statusSend[2];
  ierr = MPI_Waitall( numNei, reqRecv, statusRecv); CHKERRQ(ierr);
  ierr = MPI_Waitall( numNei, reqSend, statusSend); CHKERRQ(ierr);
 
  int j;
  for (i = 0; i < numNei; i++) {
    ierr = PetscSynchronizedPrintf( comm, "[%d] recvbuf: ", rank ); CHKERRQ(ierr);
    int *buf = ArrayGetData( recvbufs[i] );
    for (j = 0; j < ArrayLength( recvbufs[i] ); j++) {
      ierr = PetscSynchronizedPrintf( comm, "%d ", buf[j]); CHKERRQ(ierr);
    }
    ierr = PetscSynchronizedPrintf( comm, "\n" ); CHKERRQ(ierr);
  }

  ierr = PetscSynchronizedFlush(comm); CHKERRQ(ierr);

  ierr = FiberFieldDestroy(fibers); CHKERRQ(ierr);
  ierr = DCellFinalize(); CHKERRQ(ierr);
  return 0;
}
