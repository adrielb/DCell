#include "FiberField.h"

/*#define DEBUG_ALLTOALL*/

PetscErrorCode FiberField_AddToSendbufs( FiberField field );
PetscErrorCode FiberField_Nei_Alltoall( FiberField f );
PetscErrorCode FiberField_UnpackFromRecvbufs( FiberField f );

#undef __FUNCT__
#define __FUNCT__ "FiberField_SpatiallyBalance"
PetscErrorCode FiberField_SpatiallyBalance( FiberField f )
{

  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = FiberField_AddToSendbufs( f ); CHKERRQ(ierr);
  ierr = FiberField_Nei_Alltoall( f ); CHKERRQ(ierr);
  ierr = FiberField_UnpackFromRecvbufs( f ); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "FiberField_AddToSendbufs"
PetscErrorCode FiberField_AddToSendbufs( FiberField field )
{
  int i;
  int e;
  const int len = ArrayLength(field->verts);
  const Coor lmin = field->localBounds.min;
  const Coor lmax = field->localBounds.max;
  int neiIdx;
  VertexMPI *vertmpi;
  iCoor n; // index in 3x3x3 array nei
  Vertex v;
  PetscMPIInt sendRank;
  const PetscMPIInt *neiRanks;
  PetscErrorCode ierr;

  PetscFunctionBegin;

  ierr = DMDAGetNeighbors(field->da, &neiRanks); CHKERRQ(ierr);
  
  // clear send arrays
  // for each vert
  //   if outside nei, err
  //   else add vert to send list

  for (i = 0; i < NUMNEI; i++) {
    ArraySetSize( field->sendbufs[i], 0);
  }
   
  for (i = 0; i < len; i++) {
    ierr = ArrayGet( field->verts, i, &v ); CHKERRQ(ierr);

    // convert vertex Coor (v->X) to nei iCoor (n)
    // n = floor( (X-min) / (max-min) ) + 1
    // n = (1,1,1) is this local proc in 3x3x3 
    n.x = floor((v->X.x - lmin.x) / (lmax.x-lmin.x)) + 1;
    n.y = floor((v->X.y - lmin.y) / (lmax.y-lmin.y)) + 1;
    n.z = floor((v->X.z - lmin.z) / (lmax.z-lmin.z)) + 1;

    // if vertex outside 3x3x3 nei, something went terribly wrong
    if (n.x < 0 || n.x > 2 ||
        n.y < 0 || n.y > 2 ||
        n.z < 0 || n.z > 2 ) {
      ierr = PetscInfo(0, "ERROR: Vertex outside 3x3x3 neighbor region\n"); CHKERRQ(ierr);
      ierr = PetscInfo3(0, "X = {%f, %f, %f}\n",v->X.x,v->X.y,v->X.z); CHKERRQ(ierr);
      ierr = PetscInfo3(0, "n = {%d, %d, %d}\n",n.x,n.y,n.z); CHKERRQ(ierr);
      ierr = PetscInfo(0, "ERROR: END MESSAGE\n"); CHKERRQ(ierr);
      SETERRQ(field->comm, 0, "Vertex outside 3x3x3 neighbor region");
    } else { 
      // convert 3D nei coor to 1D nei index
      neiIdx = n.x + 3*n.y + 9*n.z;
      // convert nei index to mpi rank
      sendRank = neiRanks[neiIdx];
      // in the edge case where a vertex leaves the global bounding box, abort
      // handle this case in the physics, not in the communication routine
      if ( sendRank == MPI_PROC_NULL) {
        ierr = PetscInfo(0, "ERROR: Vertex outside global bbox\n"); CHKERRQ(ierr);
        ierr = PetscInfo3(0, "X = {%f, %f, %f}\n",v->X.x,v->X.y,v->X.z); CHKERRQ(ierr);
        ierr = PetscInfo3(0, "n = {%d, %d, %d}\n",n.x,n.y,n.z); CHKERRQ(ierr);
        ierr = PetscInfo1(0, "neiIdx = %d\n",neiIdx); CHKERRQ(ierr);
        ierr = PetscInfo(0, "ERROR: END MESSAGE\n"); CHKERRQ(ierr);
        SETERRQ(field->comm, 0, "Vertex outside global bbox\n");
      }

      // add vertex to send list[rank]
      ierr = ArrayAppend( field->sendbufs[neiIdx], &vertmpi); CHKERRQ(ierr); 
      vertmpi->vID = v->vID;
      vertmpi->type= v->type;
      vertmpi->X = v->X;
      vertmpi->V = v->V;
      for (e = 0; e < MAXEDGES; e++) {
        vertmpi->eID[e] = v->eID[e];

        eLO = v->eLO[e];
        edgeRank[ eLO ] = sendRank;
      }
    }
  }

  for (i = 0; i < elen; i++) {
    if (e.v0 == e.v1) {
      ierr = ArrayAppend( field->sendbufs[neiIdx], &vertmpi); CHKERRQ(ierr); 
    } else {
      ierr = ArrayAppend( field->sendbufs[neiIdx], &vertmpi); CHKERRQ(ierr); 
      ierr = ArrayAppend( field->sendbufs[neiIdx], &vertmpi); CHKERRQ(ierr); 
    }
  }

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "FiberField_Nei_Alltoall"
PetscErrorCode FiberField_Nei_Alltoall( FiberField f )
{
  int i;
  int neiIdx; // index where nei[] == src
  int count; // probing count num elements received
  const int tag = 128456826; // TODO: should tag # be something unique for each call to this routine?
  const int NUMRECV = f->NUMRECV;
  Array *sendbufs = f->sendbufs;
  Array *recvbufs = f->recvbufs;
  MPI_Request reqSend[NUMNEI];
  MPI_Request reqRecv[NUMNEI];
  MPI_Status status;
  MPI_Comm comm = f->comm;
  const PetscMPIInt *neiRanks;
  PetscErrorCode ierr;

  PetscFunctionBegin;

  //TODO: why is this barrier necessary if all Isend/Irecv matched with WaitAll?
  //BUG: without barrier, sources from other iterations caught in probe
  ierr = PetscBarrier(0); CHKERRQ(ierr);

  ierr = DMDAGetNeighbors(f->da, &neiRanks); CHKERRQ(ierr);

  // send verts to neighbors
  for (i = 0; i < NUMNEI; i++) {
    count = ArrayLength(sendbufs[i]);
    ierr = MPI_Isend(ArrayGetData(sendbufs[i]), count, f->vertmpitype, neiRanks[i], tag, comm, &reqSend[i] ); CHKERRQ(ierr);
    /*ierr = MPI_Send(ArrayGetData(sendbufs[i]), count, f->vertmpitype, neiRanks[i], tag, comm ); CHKERRQ(ierr);*/

#ifdef DEBUG_ALLTOALL
    ierr = PetscInfo1(0, "i = %d\n", i ); CHKERRQ(ierr);
    ierr = PetscInfo1(0, "dst = %d\n", neiRanks[i] ); CHKERRQ(ierr);
    ierr = PetscInfo1(0, "count = %d\n", ArrayLength(sendbufs[i]) ); CHKERRQ(ierr);
#endif
  }
  
  // receive verts from neighbors
  for (i = 0; i < NUMRECV; i++) {
    // probe for count verts sent
    ierr = MPI_Probe(MPI_ANY_SOURCE, tag, comm, &status); CHKERRQ(ierr);
    ierr = MPI_Get_count( &status, f->vertmpitype, &count); CHKERRQ(ierr);
    // convert source rank into nei index (for recvbufs array)
    for (neiIdx = 0; neiIdx < NUMNEI; neiIdx++) {
      if( neiRanks[neiIdx] == status.MPI_SOURCE )
        break;
    }

#ifdef DEBUG_ALLTOALL
    ierr = PetscInfo1(0, "i = %d\n", i ); CHKERRQ(ierr);
    ierr = PetscInfo1(0, "src = %d\n", status.MPI_SOURCE ); CHKERRQ(ierr);
    ierr = PetscInfo1(0, "count = %d\n", count ); CHKERRQ(ierr);
#endif
    
    ierr = ArraySetSize( recvbufs[neiIdx], count); CHKERRQ(ierr);
    ierr = MPI_Irecv( ArrayGetData(recvbufs[neiIdx]), count, f->vertmpitype, status.MPI_SOURCE, tag, comm, &reqRecv[i]); CHKERRQ(ierr);
    /*ierr = MPI_Recv( ArrayGetData(recvbufs[neiIdx]), count, f->vertmpitype, status.MPI_SOURCE, tag, comm, &status ); CHKERRQ(ierr);*/
  }

  ierr = MPI_Waitall( NUMRECV, reqRecv, MPI_STATUSES_IGNORE ); CHKERRQ(ierr);
  ierr = MPI_Waitall(  NUMNEI, reqSend, MPI_STATUSES_IGNORE ); CHKERRQ(ierr);
  

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "FiberField_UnpackFromRecvbufs"
PetscErrorCode FiberField_UnpackFromRecvbufs( FiberField f )
{
  int i;
  int e;
  int b;
  int l;
  int len;
  int total_recv;
  Array *recvbufs = f->recvbufs;
  Vertex v_local;
  VertexMPI *v_buf;
  PetscErrorCode ierr;

  PetscFunctionBegin;
 
  // sum up total number of verts received
  total_recv = 0;
  for (b = 0; b < NUMNEI; b++) {
    total_recv += ArrayLength( recvbufs[b] );
  }

  ierr = PetscInfo1(0, "total_recv = %d\n", total_recv  ); CHKERRQ(ierr);

  // set local vert array length to sum of total verts received
  ierr = ArraySetSize( f->verts, total_recv); CHKERRQ(ierr);
  v_local = ArrayGetData( f->verts );

  // copy from bufs to local vert list
  // for each buf
  //   for each elem in buf
  //     v_local = v_buf
  l = 0;
  for (b = 0; b < NUMNEI; b++) {
    len = ArrayLength( recvbufs[b] );
    v_buf = ArrayGetData( recvbufs[b] );
    for (i = 0; i < len; i++) {
      v_local[l].vID = v_buf[i].vID;
      v_local[l].type= v_buf[i].type;
      v_local[l].X   = v_buf[i].X;
      v_local[l].V   = v_buf[i].V;
      for (e = 0; e < MAXEDGES; e++) {
        v_local.eID[e] = v_buf.eID[e];
      }
      l++;
    }
  }

  PetscFunctionReturn(0);
}

#undef DEBUG_ALLTOALL
