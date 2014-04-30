#include "FiberField.h"

#undef __FUNCT__
#define __FUNCT__ "VertexCreate"
PetscErrorCode VertexCreate(FiberField field, Vertex *vertex)
{
  Vertex v=0;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = MemCacheAlloc(field->mcVerticies,&v); CHKERRQ(ierr);
  ierr = ArrayAppendPtr(field->fibers,v); CHKERRQ(ierr);
  ierr = UniqueIDGenerate(field->vid,&v->id); CHKERRQ(ierr);

  *vertex = v;
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "VertexDestroy"
PetscErrorCode VertexDestroy(FiberField field, Vertex vertex)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = MemCacheFree(field->mcVerticies, vertex); CHKERRQ(ierr);
  //set edge ids to zero? neg?
  //delete the edges from vertices attached to it
  //delete the edges of this vertex
  //remove from field array
  //Delete the vertex from the cache
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "VertexAddEdge"
PetscErrorCode VertexAddEdge( Vertex v0, Vertex v1, EdgeType etype )
{
  int i;

  PetscFunctionBegin;
  for (i = 0; i < MAXEDGES; ++i) {
    if( v0->v[i] == NULL ) {
      v0->v[i] = v1;
      v0->vID[i] = v1->id;
      v0->e[i] = etype;
      break;
    }
  }
  if( i == MAXEDGES ) {
    printf("MAXEDGES (%d) reached\n", MAXEDGES);
    exit(1);
  }
  PetscFunctionReturn(0);
}
