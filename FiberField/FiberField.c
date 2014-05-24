#include "FiberField.h"

#undef __FUNCT__
#define __FUNCT__ "FiberFieldCreate"
PetscErrorCode FiberFieldCreate(MPI_Comm comm, FiberField *fibers)
{
  FiberField f;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscNew(struct _FiberField, &f); CHKERRQ(ierr);
  ierr = ArrayCreate("verts",sizeof(struct _Vertex),&f->verts); CHKERRQ(ierr);
  ierr = ArrayCreate("edges",sizeof(struct _Vertex),&f->edges); CHKERRQ(ierr);
  ierr = UniqueIDCreate( &f->vid ); CHKERRQ(ierr);
  ierr = FiberField_CreateVertexMPIDatatype( f ); CHKERRQ(ierr);



  f->comm = comm;
  f->drag = 1;
  f->mass = 1;

  int i;
#define LEN 12
  char label[LEN];
  for (i = 0; i < NUMNEI; i++) {
    ierr = PetscSNPrintf( label, LEN, "sendbuf%d", i); CHKERRQ(ierr);
    ierr = ArrayCreate( label, sizeof(struct _VertexMPI), &f->sendbufs[i] ); CHKERRQ(ierr);

    ierr = PetscSNPrintf( label, LEN, "recvbuf%d", i); CHKERRQ(ierr);
    ierr = ArrayCreate( label, sizeof(struct _VertexMPI), &f->recvbufs[i] ); CHKERRQ(ierr);
  }
#undef LEN

  *fibers = f;
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "FiberFieldDestroy"
PetscErrorCode FiberFieldDestroy(FiberField fibers)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = MPI_Type_free( &fibers->vertmpitype ); CHKERRQ(ierr);
  ierr = UniqueIDDestroy(fibers->vid); CHKERRQ(ierr);
  ierr = ArrayDestroy(fibers->verts); CHKERRQ(ierr);
  ierr = ArrayDestroy(fibers->edges); CHKERRQ(ierr);
  ierr = PetscFree( fibers ); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

// FiberFieldSetup
// Use a DA to determine spatial partitioning
// TODO: DA needs to come from FluidField, not here.
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

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "FiberFieldPrint"
PetscErrorCode FiberFieldPrint( FiberField fibers )
{
  int i,j;
  int len = ArrayLength( fibers->verts );
  PetscErrorCode ierr;

  PetscFunctionBegin;
  struct _Vertex* vs = ArrayGetData(fibers->verts);
  for (i = 0; i < len; i++) {
    ierr = PetscSynchronizedPrintf( fibers->comm, "ID: %d; ", vs[i].vID ); CHKERRQ(ierr);
    for (j = 0; j < MAXEDGES; j++) {
      /*printf("%d ", vs[i].vID[j] ); */
      ierr = PetscSynchronizedPrintf( fibers->comm, "%d ", vs[i].eID[j] ); CHKERRQ(ierr);
    }
    ierr = PetscSynchronizedPrintf( fibers->comm, "\n" ); CHKERRQ(ierr);
  }
  ierr = PetscSynchronizedFlush(fibers->comm); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "FiberField_CreateVertexMPIDatatype"
PetscErrorCode FiberField_CreateVertexMPIDatatype( FiberField f )
{
  /*PetscErrorCode ierr;*/

  PetscFunctionBegin;

#define LEN 4
  MPI_Aint disp[LEN] = {
    offsetof(struct _VertexMPI, vID), // vID
    offsetof(struct _VertexMPI, eID), // eID
    offsetof(struct _VertexMPI, X),   // Coor x
    offsetof(struct _VertexMPI, V)    // Coor v
  };

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
  
  MPI_Type_create_struct(LEN, blocklen, disp, types, &f->vertmpitype );
#undef LEN

  MPI_Type_commit( &f->vertmpitype );

  PetscFunctionReturn(0);
}
