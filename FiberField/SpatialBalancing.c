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

void PositionToNeiIdx( const BoundingBox *bbox, const Coor *X, iCoor *n, int *neiIdx )
{
    // convert vertex Coor (v->X) to nei iCoor (n)
    // n = floor( (X-min) / (max-min) ) + 1
    // n == (1,1,1) is this local proc in 3x3x3 
    n->x = floor((X->x - bbox->min.x) / (bbox->max.x - bbox->min.x)) + 1;
    n->y = floor((X->y - bbox->min.y) / (bbox->max.y - bbox->min.y)) + 1;
    n->z = floor((X->z - bbox->min.z) / (bbox->max.z - bbox->min.z)) + 1;

    // convert 3D nei coor to 1D nei index
    *neiIdx = n->x + 3*n->y + 9*n->z;
}

#undef __FUNCT__
#define __FUNCT__ "FiberField_AddToSendbufs"
PetscErrorCode FiberField_AddToSendbufs( FiberField field )
{
  int i;
  int e;
  const int vlen = ArrayLength(field->verts);
  const int elen = ArrayLength(field->edges);
  const BoundingBox lbbox = field->localBounds;
  int neiIdx;
  VertexEdgeMPI *evmpi;
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
   
  for (i = 0; i < vlen; i++) {
    ierr = ArrayGet( field->verts, i, &v ); CHKERRQ(ierr);

    PositionToNeiIdx( &lbbox, &v->X, &n, &neiIdx);

    // if vertex outside 3x3x3 nei, something went terribly wrong
    if (n.x < 0 || n.x > 2 ||
        n.y < 0 || n.y > 2 ||
        n.z < 0 || n.z > 2 ) {
      ierr = PetscInfo(0, "ERROR: Vertex outside 3x3x3 neighbor region\n"); CHKERRQ(ierr);
      ierr = PetscInfo1(0, "i = %d\n",i); CHKERRQ(ierr);
      ierr = PetscInfo3(0, "X = {%f, %f, %f}\n",v->X.x,v->X.y,v->X.z); CHKERRQ(ierr);
      ierr = PetscInfo3(0, "n = {%d, %d, %d}\n",n.x,n.y,n.z); CHKERRQ(ierr);
      ierr = PetscInfo(0, "ERROR: END MESSAGE\n"); CHKERRQ(ierr);
      SETERRQ(field->comm, 0, "Vertex outside 3x3x3 neighbor region");
    } else { 
      // convert nei index to mpi rank
      sendRank = neiRanks[neiIdx];
      // in the edge case where a vertex leaves the global bounding box, abort
      // handle this case in the physics, not in the communication routine
      if ( sendRank == MPI_PROC_NULL) {
        ierr = PetscInfo(0, "ERROR: Vertex outside global bbox\n"); CHKERRQ(ierr);
        ierr = PetscInfo1(0, "i = %d\n",i); CHKERRQ(ierr);
        ierr = PetscInfo3(0, "X = {%f, %f, %f}\n",v->X.x,v->X.y,v->X.z); CHKERRQ(ierr);
        ierr = PetscInfo3(0, "n = {%d, %d, %d}\n",n.x,n.y,n.z); CHKERRQ(ierr);
        ierr = PetscInfo1(0, "neiIdx = %d\n",neiIdx); CHKERRQ(ierr);
        ierr = PetscInfo(0, "ERROR: END MESSAGE\n"); CHKERRQ(ierr);
        SETERRQ(field->comm, 0, "Vertex outside global bbox\n");
      }

      // add vertex to send list[rank]
      ierr = ArrayAppend( field->sendbufs[neiIdx], &evmpi); CHKERRQ(ierr); 
      evmpi->xID = v->vID;
      evmpi->type= v->type;
      evmpi->X = v->X;
      evmpi->V = v->V;
      for (e = 0; e < MAXEDGES; e++) {
        evmpi->yIDs[e] = v->eID[e];
      }
    }
  }

  int min;
  int vPO;
  struct _Edge *edges = ArrayGetData(field->edges);
  struct _Vertex *vertsPO;
  ierr = FiberFieldGetVertexArrayPO( field, &vertsPO ); CHKERRQ(ierr);
  for (e = 0; e < elen; e++) {
    // the edge is 'owned' by the vertex with the smallest ID
    min = edges[e].vID[0] < edges[e].vID[1] ? 0 : 1;
    vPO = edges[e].vPO[min]; 

    v = &vertsPO[vPO];

    PositionToNeiIdx( &lbbox, &v->X, &n, &neiIdx);

    if (v->vID != edges[e].vID[min] ) {
      ierr = PetscInfo1(0, "v->vID = %d\n", v->vID); CHKERRQ(ierr);
      ierr = PetscInfo1(0, "edges[e].vID[min] = %d\n", edges[e].vID[min]); CHKERRQ(ierr);
      SETERRQ(PETSC_COMM_SELF, 0, "Bad vertex");
    }

    ierr = ArrayAppend( field->sendbufs[neiIdx], &evmpi); CHKERRQ(ierr); 
    evmpi->xID = edges[e].eID; 
    evmpi->type = edges[e].type; 
    evmpi->yIDs[0] = edges[e].vID[0]; 
    evmpi->yIDs[1] = edges[e].vID[1]; 
    evmpi->X.x = edges[e].l0;
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
  int lv; // local vert index
  int le; // local edge index
  int len; // length of recvbuf
  int total_recv_verts;
  int total_recv_edges;
  Array *recvbufs = f->recvbufs;
  Vertex v_local;
  Edge e_local;
  VertexEdgeMPI *buf;
  FiberType *ftype = ArrayGetData( f->fibertypesDB );
  PetscErrorCode ierr;

  PetscFunctionBegin;
 
  // sum up total number of verts and edges received
  total_recv_verts = 0;
  total_recv_edges = 0;

  for (b = 0; b < NUMNEI; b++) {
    len = ArrayLength( recvbufs[b] );
    buf = ArrayGetData( recvbufs[b] );
    for (i = 0; i < len; i++) {
      if ( ftype[ buf[i].type ].isEdge ) {
        total_recv_edges++;
      } else {
        total_recv_verts++;
      }
    }
  }

  ierr = PetscInfo1(0, "total_recv_verts = %d\n", total_recv_verts  ); CHKERRQ(ierr);
  ierr = PetscInfo1(0, "total_recv_edges = %d\n", total_recv_edges  ); CHKERRQ(ierr);

  // set local vert/edge array length to sum of total received
  ierr = ArraySetSize( f->verts, total_recv_verts); CHKERRQ(ierr);
  ierr = ArraySetSize( f->edges, total_recv_edges); CHKERRQ(ierr);
  v_local = ArrayGetData( f->verts );
  e_local = ArrayGetData( f->edges );

  // copy from bufs to local vert list
  // for each buf
  //   for each elem in buf
  //     if edge
  //       e_local = buf
  //     else
  //       v_local = buf
  lv = 0;
  le = 0;
  for (b = 0; b < NUMNEI; b++) {
    len = ArrayLength( recvbufs[b] );
    buf = ArrayGetData( recvbufs[b] );
    for (i = 0; i < len; i++) {
      if ( ftype[ buf[i].type ].isEdge ) {
        e_local[le].eID    = buf[i].xID;
        e_local[le].type   = buf[i].type;
        e_local[le].l0     = buf[i].X.x;
        e_local[le].vID[0] = buf[i].yIDs[0];
        e_local[le].vID[1] = buf[i].yIDs[1];
        le++;
      } else {
        v_local[lv].vID = buf[i].xID;
        v_local[lv].type= buf[i].type;
        v_local[lv].X   = buf[i].X;
        v_local[lv].V   = buf[i].V;
        for (e = 0; e < MAXEDGES; e++) {
          v_local[lv].eID[e] = buf[i].yIDs[e];
        }
        lv++;
      } // if vertex or edge
    } // for i in recvbuf
  } // for b in NUMNEI

  PetscFunctionReturn(0);
}

#undef DEBUG_ALLTOALL
