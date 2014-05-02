#ifndef FIBERFIELD_H_
#define FIBERFIELD_H_

#include "Common.h"

typedef struct _FiberField *FiberField;
typedef struct _Vertex *Vertex;
typedef int VertexType;
typedef int EdgeType;
typedef UniqueIDType VertexID;

#define MAXEDGES 4
#define FIBERFIELD_NO_VERTEX 0

struct _Vertex {
  VertexID ID;     // Unique ID for this vertex
  VertexType type; // type of vertex
  Coor X;          // position
  Coor V;          // velocity
  // Connected vertices 
  Vertex      v[MAXEDGES];
  VertexID  vID[MAXEDGES];
  EdgeType    e[MAXEDGES];
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
PetscErrorCode FiberFieldPrint( FiberField fibers );
PetscErrorCode VertexCreate(FiberField field, Vertex *v);
PetscErrorCode VertexAddEdge( Vertex v0, Vertex v1, EdgeType etype );
PetscErrorCode VertexRemoveEdge( Vertex v0, Vertex v1 );

PetscErrorCode Vertex_Link( Vertex v0, Vertex v1, EdgeType etype );
PetscErrorCode Vertex_Unlink( Vertex v0, VertexID v1_id );
#endif /* FIBERFIELD_H_ */
