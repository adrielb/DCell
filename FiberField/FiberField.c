#include "FiberField.h"
#include "FiberField_private.h"

PetscErrorCode FiberField_CreateVertexMPIDatatype( FiberField f );

#undef __FUNCT__
#define __FUNCT__ "FiberFieldCreate"
PetscErrorCode FiberFieldCreate(MPI_Comm comm, FiberField *fibers)
{
  FiberField f;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscNew(struct _FiberField, &f); CHKERRQ(ierr);
  ierr = ArrayCreate("verts",sizeof(struct _Vertex),&f->verts); CHKERRQ(ierr);
  ierr = ArrayCreate("edges",sizeof(struct _Edge),&f->edges); CHKERRQ(ierr);
  ierr = ArrayCreate("vIDs",sizeof(VertexID),&f->vIDs); CHKERRQ(ierr);
  ierr = ArrayCreate("eIDs",sizeof(EdgeID),  &f->eIDs); CHKERRQ(ierr);
  ierr = ArrayCreate("fiberDB", sizeof(FiberType), &f->fibertypesDB); CHKERRQ(ierr);
  ierr = UniqueIDCreate( &f->vid ); CHKERRQ(ierr);
  ierr = UniqueIDCreate( &f->eid ); CHKERRQ(ierr);
  ierr = FiberField_CreateVertexMPIDatatype( f ); CHKERRQ(ierr);
  ierr = PetscRandomCreate(PETSC_COMM_WORLD,&f->rnd); CHKERRQ(ierr);
  ierr = PetscRandomSetType(f->rnd,PETSCRAND48); CHKERRQ(ierr);
  /*ierr = PetscRandomSetFromOptions(rnd); CHKERRQ(ierr);*/
  /*ierr = PetscRandomSeed(rnd); CHKERRQ(ierr);*/
  ierr = ArrayCreate("fiber", sizeof(Vertex), &f->fiber); CHKERRQ(ierr);
  

  ierr = FiberFieldSetFluidVelocityEvaluator( f, FiberField_ZeroFluidVelocity ); CHKERRQ(ierr);


  f->comm = comm;
  f->fluidDrag = 1;
  f->mass = 1;
  f->dt = 0;

  int i;
  char label[13];
  for (i = 0; i < NUMNEI; i++) {
    ierr = PetscSNPrintf( label, sizeof(label), "sendbuf%d", i); CHKERRQ(ierr);
    ierr = ArrayCreate( label, sizeof(struct _VertexEdgeMPI), &f->sendbufs[i] ); CHKERRQ(ierr);

    ierr = PetscSNPrintf( label, sizeof(label), "recvbuf%d", i); CHKERRQ(ierr);
    ierr = ArrayCreate( label, sizeof(struct _VertexEdgeMPI), &f->recvbufs[i] ); CHKERRQ(ierr);
  }

  *fibers = f;
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "FiberFieldDestroy"
PetscErrorCode FiberFieldDestroy(FiberField fibers)
{
  int i;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscInfo(0, "Destroying Fiber Field\n"); CHKERRQ(ierr);
  ierr = ArrayDestroy(fibers->verts); CHKERRQ(ierr);
  ierr = ArrayDestroy(fibers->edges); CHKERRQ(ierr);
  ierr = ArrayDestroy(fibers->vIDs); CHKERRQ(ierr);
  ierr = ArrayDestroy(fibers->eIDs); CHKERRQ(ierr);
  ierr = UniqueIDDestroy(fibers->vid); CHKERRQ(ierr);
  ierr = UniqueIDDestroy(fibers->eid); CHKERRQ(ierr);
  ierr = MPI_Type_free( &fibers->vertmpitype ); CHKERRQ(ierr);
  ierr = PetscRandomDestroy(&fibers->rnd); CHKERRQ(ierr);
  ierr = ArrayDestroy(fibers->fiber); CHKERRQ(ierr);
  for (i = 0; i < NUMNEI; i++) {
    ierr = ArrayDestroy( fibers->sendbufs[i] ); CHKERRQ(ierr);
    ierr = ArrayDestroy( fibers->recvbufs[i] ); CHKERRQ(ierr);
  }
  ierr = PetscFree( fibers ); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

// FiberFieldSetup
// Use a DA to determine spatial partitioning
// TODO: DA needs to come from FluidField, not here.
// TODO: DA needs to match MPI_Cart_create
#undef __FUNCT__
#define __FUNCT__ "FiberFieldSetup"
PetscErrorCode FiberFieldSetup( FiberField fibers )
{
  DMDALocalInfo info;
  const PetscReal dh = fibers->dh;
  const int dof = 1;
  const Coor gmin = fibers->globalBounds.min;
  const Coor gmax = fibers->globalBounds.max;
  Coor *dX = &fibers->dX;
  PetscErrorCode ierr;

  PetscFunctionBegin;

  dX->x = dh; 
  dX->y = dh; 
  dX->z = dh; 

  iCoor dims = { 
    (gmax.x - gmin.x) / dX->x,
    (gmax.y - gmin.y) / dX->y, 
    (gmax.z - gmin.z) / dX->z
  };

  ierr = DMDACreate3d(fibers->comm,
      DMDA_BOUNDARY_NONE,DMDA_BOUNDARY_NONE,DMDA_BOUNDARY_NONE,
      DMDA_STENCIL_BOX,
      dims.x,dims.y,dims.z, 
      PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE,
      dof,1,  0,0,0, &fibers->da); CHKERRQ(ierr);

  ierr = DMDAGetLocalInfo(fibers->da, &info); CHKERRQ(ierr);

  BoundingBox localBounds = {
    .min = {
      gmin.x + dX->x * info.xs,
      gmin.y + dX->y * info.ys,
      gmin.z + dX->z * info.zs
    }, 
    .max = {
      gmin.x + dX->x * (info.xs + info.xm),
      gmin.y + dX->y * (info.ys + info.ym),
      gmin.z + dX->z * (info.zs + info.zm)
    }
  };
  fibers->localBounds = localBounds;

  // Do not wait for receives from NULL rank neighbors
  int i;
  const int *neiRanks;
  ierr = DMDAGetNeighbors(fibers->da, &neiRanks); CHKERRQ(ierr);
  fibers->NUMRECV = NUMNEI;
  char strNei[64] = {0};
  char strRank[4];
  for (i = 0; i < NUMNEI; ++i) {
    if (neiRanks[i] == MPI_PROC_NULL) {
      fibers->NUMRECV--;
    } else {
      sprintf( strRank, " %d", neiRanks[i] );
      strcat( strNei, strRank );
    }
  }

  ierr = PetscInfo3(fibers->da, "globalBounds.min = %f, %f, %f\n", gmin.x, gmin.y, gmin.z ); CHKERRQ(ierr);
  ierr = PetscInfo3(fibers->da, "globalBounds.max = %f, %f, %f\n", gmax.x, gmax.y, gmax.z ); CHKERRQ(ierr);
  ierr = PetscInfo3(0, "localBounds.min = %f, %f, %f\n", localBounds.min.x, localBounds.min.y, localBounds.min.z ); CHKERRQ(ierr);
  ierr = PetscInfo3(0, "localBounds.max = %f, %f, %f\n", localBounds.max.x, localBounds.max.y, localBounds.max.z ); CHKERRQ(ierr);
  ierr = PetscInfo1(0, "NUMRECV = %d\n", fibers->NUMRECV ); CHKERRQ(ierr);
  ierr = PetscInfo1(0, "strNei =%s\n", strNei ); CHKERRQ(ierr);

  int len = ArrayLength( fibers->fibertypesDB );
  FiberType *ftypes = ArrayGetData( fibers->fibertypesDB );
  for (i = 0; i < len; i++) {
    ierr = PetscInfo1(fibers->da, "ID = %d\n", ftypes[i].ID); CHKERRQ(ierr);
    ierr = PetscInfo1(fibers->da, "isEdge = %d\n", ftypes[i].isEdge); CHKERRQ(ierr);
    ierr = PetscInfo1(fibers->da, "name = %s\n", ftypes[i].name); CHKERRQ(ierr);
  }
  if (len == 0) {
    SETERRQ(PETSC_COMM_WORLD,0,"FiberTypes DB not set\n");
  }

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "FiberFieldPrint"
PetscErrorCode FiberFieldPrint( FiberField fibers )
{
  int i,j;
  const int vlen = ArrayLength( fibers->verts );
  const int elen = ArrayLength( fibers->edges );
  const struct _Vertex* vs = ArrayGetData(fibers->verts);
  const struct _Edge* edges = ArrayGetData(fibers->edges);
  PetscErrorCode ierr;

  PetscFunctionBegin;
  for (i = 0; i < vlen; i++) {
    ierr = PetscSynchronizedPrintf( fibers->comm, "vID: %d; ", vs[i].vID ); CHKERRQ(ierr);
    for (j = 0; j < MAXEDGES; j++) {
      ierr = PetscSynchronizedPrintf( fibers->comm, "%d ", vs[i].eID[j] ); CHKERRQ(ierr);
    }
    ierr = PetscSynchronizedPrintf( fibers->comm, "\n" ); CHKERRQ(ierr);
  }
  ierr = PetscSynchronizedFlush(fibers->comm); CHKERRQ(ierr);

  for (i = 0; i < elen; i++) {
    ierr = PetscSynchronizedPrintf( fibers->comm, "eID: %d; %d %d\n", edges[i].eID, edges[i].vID[0], edges[i].vID[1] ); CHKERRQ(ierr);
  }
  ierr = PetscSynchronizedFlush(fibers->comm); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "FiberFieldView"
PetscErrorCode FiberFieldView( FiberField f )
{
  int i;
  int len;
  PetscReal *x0, *xi, *x1, *Fe_v;
  PetscReal *v0, *vi, *v1;
  MPI_Comm comm = f->comm;
  PetscErrorCode ierr;

  PetscFunctionBegin;

  ierr = PetscPrintf(PETSC_COMM_WORLD, "\nAO Verts\n"); CHKERRQ(ierr);
  ierr = AOView( f->aoVerts, PETSC_VIEWER_STDOUT_WORLD ); CHKERRQ(ierr);

  ierr = PetscPrintf(PETSC_COMM_WORLD, "\nAO Edges\n"); CHKERRQ(ierr);
  ierr = AOView( f->aoEdges, PETSC_VIEWER_STDOUT_WORLD ); CHKERRQ(ierr);

  ierr = PetscPrintf(PETSC_COMM_WORLD, "\nX to U\n"); CHKERRQ(ierr);
  ierr = MatView( f->matXtoU, PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);

  ierr = PetscPrintf(PETSC_COMM_WORLD, "\nFe to Fv\n"); CHKERRQ(ierr);
  ierr = MatView( f->matFe2Fv, PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);

  ierr = VecGetLocalSize(f->x0, &len); CHKERRQ(ierr);
  ierr = VecGetArray( f->x0, &x0); CHKERRQ(ierr);
  ierr = VecGetArray( f->xi, &xi); CHKERRQ(ierr);
  ierr = VecGetArray( f->x1, &x1); CHKERRQ(ierr);
  ierr = VecGetArray( f->v0, &v0); CHKERRQ(ierr);
  ierr = VecGetArray( f->vi, &vi); CHKERRQ(ierr);
  ierr = VecGetArray( f->v1, &v1); CHKERRQ(ierr);
  ierr = VecGetArray( f->Fe_v, &Fe_v); CHKERRQ(ierr);
  ierr = PetscPrintf(comm, "x0,v0,xi,vi,x1,v1,Fe_v\n" ); CHKERRQ(ierr);
  for (i = 0; i < len; i++) {
    ierr = PetscSynchronizedPrintf(comm, "%f,%f,%f,%f,%f,%f,%f\n", x0[i], v0[i], 
                                                                   xi[i], vi[i],
                                                                   x1[i], v1[i],
                                                                   Fe_v[i] ); CHKERRQ(ierr);
  }
  ierr = PetscSynchronizedFlush(comm); CHKERRQ(ierr);
  ierr = VecRestoreArray( f->x0, &x0); CHKERRQ(ierr);
  ierr = VecRestoreArray( f->xi, &xi); CHKERRQ(ierr);
  ierr = VecRestoreArray( f->x1, &x1); CHKERRQ(ierr);
  ierr = VecRestoreArray( f->v0, &v0); CHKERRQ(ierr);
  ierr = VecRestoreArray( f->vi, &vi); CHKERRQ(ierr);
  ierr = VecRestoreArray( f->v1, &v1); CHKERRQ(ierr);
  ierr = VecRestoreArray( f->Fe_v, &Fe_v); CHKERRQ(ierr);

  /*ierr = PetscPrintf(PETSC_COMM_WORLD, "\nx0\n"); CHKERRQ(ierr);*/
  /*ierr = VecView( f->x0, PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);*/

  /*ierr = PetscPrintf(PETSC_COMM_WORLD, "\nxi\n"); CHKERRQ(ierr);*/
  /*ierr = VecView( f->xi, PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);*/

  /*ierr = PetscPrintf(PETSC_COMM_WORLD, "\nx1\n"); CHKERRQ(ierr);*/
  /*ierr = VecView( f->x1, PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);*/

  /*ierr = PetscPrintf(PETSC_COMM_WORLD, "\nFe_v\n"); CHKERRQ(ierr);*/
  /*ierr = VecView( f->Fe_v, PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);*/

  ierr = PetscPrintf(PETSC_COMM_WORLD, "\nFe_e\n"); CHKERRQ(ierr);
  ierr = VecView( f->Fe_e, PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);

  /*ierr = PetscPrintf(PETSC_COMM_WORLD, "\nv0\n"); CHKERRQ(ierr);*/
  /*ierr = VecView( f->v0, PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);*/

  /*ierr = PetscPrintf(PETSC_COMM_WORLD, "\nvi\n"); CHKERRQ(ierr);*/
  /*ierr = VecView( f->vi, PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);*/

  /*ierr = PetscPrintf(PETSC_COMM_WORLD, "\nv1\n"); CHKERRQ(ierr);*/
  /*ierr = VecView( f->v1, PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);*/

  ierr = PetscPrintf(PETSC_COMM_WORLD, "\nu\n"); CHKERRQ(ierr);
  ierr = VecView( f->u, PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);

  ierr = PetscPrintf(PETSC_COMM_WORLD, "\nw\n"); CHKERRQ(ierr);
  ierr = VecView( f->w, PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "FiberField_CreateVertexMPIDatatype"
PetscErrorCode FiberField_CreateVertexMPIDatatype( FiberField f )
{
  /*PetscErrorCode ierr;*/

  PetscFunctionBegin;

  MPI_Aint disp[] = {
    offsetof(struct _VertexEdgeMPI, xID),  // xID
    offsetof(struct _VertexEdgeMPI, type), // type
    offsetof(struct _VertexEdgeMPI, yIDs), // yIDs
    offsetof(struct _VertexEdgeMPI, X),    // Coor x
    offsetof(struct _VertexEdgeMPI, V)     // Coor v
  };

  MPI_Datatype types[] = {
    MPI_INT,    // xID
    MPI_INT,    // type
    MPI_INT,    // yIDs
    MPI_DOUBLE, // Coor x
    MPI_DOUBLE  // Coor v
  };

  int blocklen[] = { 
    1,        // xID
    1,        // type
    MAXEDGES, // yIDs
    3,        // Coor x
    3         // Coor v
  };

  int count = sizeof(blocklen) / sizeof(int);

  MPI_Type_create_struct(count, blocklen, disp, types, &f->vertmpitype );

  MPI_Type_commit( &f->vertmpitype );

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "FiberFieldWrite"
PetscErrorCode FiberFieldWrite( FiberField f, int ti )
{
  int i;
  int rank; 
  char edgefilename[PETSC_MAX_PATH_LEN];
  FILE *edgefile;
  const int elen = ArrayLength( f->edges );
  struct _Edge *edges = ArrayGetData( f->edges );
  Edge e;
  PetscErrorCode ierr;

  PetscFunctionBegin;

  // handle special case of first time step,
  // values not in x1 yet until after FiberFieldSolve()
  if (ti == 0) {
    ierr = FiberField_Init( f ); CHKERRQ(ierr);
    ierr = FiberField_ToVec( f ); CHKERRQ(ierr);
    ierr = VecCopy(f->x0,f->x1); CHKERRQ(ierr);
  }

  ierr = VecWrite( f->x1, "verts", ti ); CHKERRQ(ierr);

  // Write edge list as CSV file
  MPI_Comm_rank( f->comm, &rank );
  sprintf( edgefilename, "edgelist.%d.%d.csv", ti + FILE_COUNT_START, rank);
  ierr = PetscFOpen(PETSC_COMM_SELF, edgefilename, "w", &edgefile);CHKERRQ(ierr);
  
  // Add a header line
  //if (rank == 0) {
  if( PETSC_TRUE ) {
    ierr = PetscFPrintf( PETSC_COMM_SELF, edgefile, "eID,etype,vID0,vID1,vPO0,vPO1,eP0\n"); CHKERRQ(ierr);
  }
  for (i = 0; i < elen; i++) {
    e = &edges[i];
    //                                             eI,et,vI,vI,vP,vP,eP
    ierr = PetscFPrintf(PETSC_COMM_SELF, edgefile,"%d,%d,%d,%d,%d,%d,%d\n",
        e->eID, e->type, e->vID[0], e->vID[1], e->vPO[0], e->vPO[1], e->ePO ); CHKERRQ(ierr);
  }

  ierr = PetscFClose(PETSC_COMM_SELF, edgefile); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "FiberFieldAddType"
PetscErrorCode FiberFieldAddType( FiberField f, const char *typeName, PetscBool isEdge, FiberTypeID *ID )
{
  FiberType *ftype;
  PetscErrorCode ierr;

  PetscFunctionBegin;

  *ID = ArrayLength( f->fibertypesDB );

  ierr = ArrayAppend( f->fibertypesDB, &ftype); CHKERRQ(ierr);
  ierr = PetscStrncpy( ftype->name, typeName, sizeof(ftype->name) ); CHKERRQ(ierr);
  ftype->isEdge = isEdge;
  ftype->ID = *ID;

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "FiberFieldGetVertexArrayPO"
PetscErrorCode FiberFieldGetVertexArrayPO( FiberField f, Vertex *vertsPO )
{
  int start;
  int end;
  struct _Vertex *verts;
  PetscErrorCode ierr;

  PetscFunctionBegin;

  verts = ArrayGetData( f->verts );

  /*verts[0].vPO*/

  ierr = VecGetOwnershipRange(f->x1, &start, &end); CHKERRQ(ierr);

  start = start / 3;

  *vertsPO = verts - start;

  PetscFunctionReturn(0);
}
