#ifndef FIBERFIELD_H_
#define FIBERFIELD_H_

#include "Common.h"

typedef struct _FiberField *FiberField;
typedef struct _Vertex *Vertex;
typedef struct _Edge *Edge;

#define MAXEDGES 4

struct _Vertex {
  UniqueIDType id;   // Unique ID for this vertex
  int type; // type of vertex
  Coor X; // position
  Coor V; // velocity
//  slist edgeIDs;  // list of edge IDs attached to this vertex
  UniqueIDType edgeIDs[MAXEDGES];
  Edge edges[MAXEDGES]; // list of edge pointers attached to this vertex
};

struct _Edge {
  UniqueIDType id;
  int type;
  UniqueIDType id_a;  // id of vertex A
  Vertex a;  // pointer to vertex A
  UniqueIDType id_b;  // id of vertex B
  Vertex b;  // pointer to vertex B
};

struct {
  Coor min;
  Coor max;
} BoundingBox;

struct _FiberField {
  UniqueID vid; // vertex id generator
  UniqueID eid; // edge id generator
  PetscReal mass;      // Mass of node
  PetscReal thickness; // Fiber thickness ~0.005um?
  PetscReal TOL;
  MemCache mcVerticies;
  MemCache mcEdges;
  Array fibers;
};

PetscErrorCode FiberFieldCreate(MPI_Comm comm, FiberField *fibers);
PetscErrorCode FiberFieldDestroy(FiberField fibers);
PetscErrorCode VertexCreate(FiberField field, Vertex *v);
PetscErrorCode VertexAddEdge( Vertex v, Edge e);

#endif /* FIBERFIELD_H_ */
