#ifndef FIBERFIELD_H_
#define FIBERFIELD_H_

#include "Common.h"

typedef struct _FiberField *FiberField;
typedef struct _Vertex *Vertex;
typedef int VertexType;
typedef int EdgeType;

#define MAXEDGES 4

struct _Vertex {
  UniqueIDType id; // Unique ID for this vertex
  VertexType type; // type of vertex
  Coor X;          // position
  Coor V;          // velocity
  Vertex          v[MAXEDGES];
  UniqueIDType  vID[MAXEDGES];
  EdgeType        e[MAXEDGES];
  
};

struct {
  Coor min;
  Coor max;
} BoundingBox;

struct _FiberField {
  UniqueID vid; // vertex id generator
  PetscReal mass;      // Mass of node
  PetscReal thickness; // Fiber thickness ~0.005um?
  PetscReal TOL;
  Array fibers;
};

PetscErrorCode FiberFieldCreate(MPI_Comm comm, FiberField *fibers);
PetscErrorCode FiberFieldDestroy(FiberField fibers);
PetscErrorCode VertexCreate(FiberField field, Vertex *v);
PetscErrorCode VertexAddEdge( Vertex v0, Vertex v1, EdgeType etype );
PetscErrorCode VertexRemoveEdge( Vertex v0, Vertex v1 );

#endif /* FIBERFIELD_H_ */
