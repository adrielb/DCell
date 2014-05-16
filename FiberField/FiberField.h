#ifndef FIBERFIELD_H_
#define FIBERFIELD_H_

#include "Common.h"

typedef struct _FiberField *FiberField;
typedef struct _Vertex *Vertex;
typedef struct _Edge *Edge;
typedef int VertexType;
typedef int EdgeType;
typedef UniqueIDType VertexID;
typedef UniqueIDType EdgeID;

#define MAXEDGES 4
#define FIBERFIELD_NO_VERTEX -1
#define FIBERFIELD_NO_EDGE -1

struct _Vertex {
  int petscIndex;       // Global index of edg
  VertexID vID;          // Unique ID for this vertex
  VertexType type;      // type of vertex
  Edge edges[MAXEDGES]; // Connected edges (pointers)
  EdgeID eID[MAXEDGES]; // Connected edges (ID number) 
  int    ePO[MAXEDGES]; // global petsc ordering of edge IDs
  Coor X;
};

struct _Edge {
  int petscIndex; // Global index of edge
  EdgeID eID;
  EdgeType type;
  VertexID vID[2]; // global unique ID of vertex ID
  int      vPO[2]; // global petsc ordering of vertex ID
};

struct {
  Coor min;
  Coor max;
} BoundingBox;

struct _FiberField {
  MPI_Comm comm;
  UniqueID vid;        // vertex id generator
  PetscReal mass;      // Mass of node
  PetscReal thickness; // Fiber thickness ~0.005um?
  PetscReal drag;      // drag btw fluid and vertex velocity
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
};

PetscErrorCode FiberFieldCreate(MPI_Comm comm, FiberField *fibers);
PetscErrorCode FiberFieldDestroy(FiberField fibers);
PetscErrorCode FiberFieldPrint( FiberField fibers );
PetscErrorCode VertexCreate(FiberField field, Vertex *v);
PetscErrorCode VertexAddEdge( Vertex v0, Vertex v1, EdgeType etype );
PetscErrorCode VertexRemoveEdge( Vertex v0, Vertex v1 );

PetscErrorCode Vertex_Link( Vertex v0, Vertex v1, EdgeType etype );
PetscErrorCode Vertex_Unlink( Vertex v0, VertexID v1_id );


PetscErrorCode FiberField_Setup( FiberField field );
PetscErrorCode FiberField_Step( FiberField field );
#endif /* FIBERFIELD_H_ */
