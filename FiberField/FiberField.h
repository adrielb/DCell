#ifndef FIBERFIELD_H_
#define FIBERFIELD_H_

#include "Common.h"

typedef struct _FiberField *FiberField;
typedef struct _Vertex *Vertex;
typedef struct _Edge *Edge;
typedef int FiberTypeID;
typedef UniqueIDType VertexID;
typedef UniqueIDType EdgeID;

#define NUMNEI 27
#define NEI_PROC_LOCAL 13
#define MAXEDGES 4
#define FIBERFIELD_NO_VERTEX -1
#define FIBERFIELD_NO_EDGE -1

struct _Vertex {
  VertexID vID;         // Unique ID for this vertex
  FiberTypeID type;     // type of vertex
  EdgeID eID[MAXEDGES]; // Connected edges (ID number) 
  int    ePO[MAXEDGES]; // global petsc ordering of edge IDs
  int    vPO;           // Global petsc index of vert
  Coor X;
  Coor V;
};

typedef struct _VertexEdgeMPI {
  int xID; // x=Edge, y=Verts or x=Vert, y=Edges
  FiberTypeID type; 
  int yIDs[MAXEDGES]; // IDs vertex/edge is connected too
  Coor X;
  Coor V;
} VertexEdgeMPI;

struct _Edge {
  EdgeID eID;
  FiberTypeID type;
  VertexID vID[2]; // global unique ID of vertex ID
  PetscReal l0;    // rest length of edge
  int      vPO[2]; // global petsc ordering of vertex ID
  int      ePO;    // Global index of edge
};

typedef struct {
  Coor min;
  Coor max;
} BoundingBox;

typedef struct {
  FiberTypeID ID;
  PetscBool isEdge;
  char name[32];
} FiberType;

struct _FiberField {
  MPI_Comm comm;
  MPI_Datatype vertmpitype; 
  PetscReal dh;
  Coor dX; // 
  BoundingBox globalBounds; // global bounding box of world
  BoundingBox localBounds;  // local processor bbox
  DM da;
  PetscMPIInt *neiRanks; // neiRank = DMDAGetNeighbors(da)
  int NUMRECV; // number of non-null nei ranks

  UniqueID vid;        // vertex ID generator
  UniqueID eid;        // edge ID generator
  PetscReal mass;      // Mass of node
  PetscReal thickness; // Fiber thickness ~0.005um?
  PetscReal fluidDrag;      // drag btw fluid and vertex velocity
  PetscReal TOL;
  PetscReal dt; // Time step
  PetscReal ti; // Intermediate time
  PetscReal t; // current time
  Array verts;
  Array edges;
  Mat matXtoU;     // vertex positions to displacments for each edge
  Mat matFe2Fv;    // summing edge forces for each vertex
                   // length = 3 * |verts|
  Vec x0, xi, x1;  // position at start, intermediate, end step
  Vec v0, vi, v1;  // velocity at start, intermediate, end step
  Vec Fe_v;        // aggregate elastic force on vertex
  Vec vf;          // Fluid velocity at vertex
                   // length = |edges|
  Vec l0;          // rest length
                   // length = 3 * |edges|
  Vec u;           // u_k = x_j - x_i, i < j
  Vec w;           // edge direction, unit vector
  Vec Fe_e;        // elastic force on edge, Fe_v = Sum_i( Fe_e[i] * w[i] )
  
  AO aoVerts; 
  AO aoEdges;
  Array vIDs;
  Array eIDs;

  Array sendbufs[NUMNEI];
  Array recvbufs[NUMNEI];

  Array fibertypesDB; // 

  PetscRandom rnd;
  Array fiber; // 

  PetscErrorCode (*EvaluateFluidVelocity)(FiberField f, PetscReal t, Vec x);
};

PetscErrorCode FiberFieldCreate(MPI_Comm comm, FiberField *fibers);
PetscErrorCode FiberFieldDestroy(FiberField fibers);
PetscErrorCode FiberFieldSetup(FiberField fibers);
PetscErrorCode FiberFieldPrint( FiberField fibers );
PetscErrorCode FiberFieldView( FiberField fibers );
PetscErrorCode FiberFieldWrite( FiberField fibers, int ti );

PetscErrorCode FiberFieldAddType( FiberField f, const char *typeName, PetscBool isEdge, FiberTypeID *ID );
PetscErrorCode FiberFieldGetVertexArrayPO( FiberField f, Vertex *vertsPO );
PetscErrorCode FiberFieldAddVertex(FiberField field, FiberTypeID vtype, Vertex *v);
PetscErrorCode FiberFieldAddEdge( FiberField f, Vertex v0, Vertex v1, FiberTypeID etype, PetscReal l0 );
PetscErrorCode FiberFieldRemoveEdge( Vertex v0, Vertex v1 );

PetscErrorCode FiberFieldSolve( FiberField field );

PetscErrorCode FiberFieldSetFluidVelocityEvaluator( FiberField f, void *func );
PetscErrorCode FiberFieldEvaluateFluidVelocity( FiberField f, PetscReal t, Vec x );
PetscErrorCode FiberField_ZeroFluidVelocity( FiberField f, PetscReal t, Vec x );
PetscErrorCode FiberField_CircularFluidVelocity( FiberField f, PetscReal t, Vec x );

PetscErrorCode FiberField_SpatiallyBalance( FiberField f );

PetscErrorCode FiberFieldInitLocalFiber( FiberField fibers, int numVerticies, PetscReal l0,
    FiberTypeID vertexType, FiberTypeID edgeType, FiberTypeID bendingEdgeType );

#endif /* FIBERFIELD_H_ */
