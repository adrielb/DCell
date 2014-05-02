#include "FiberField.h"

#undef __FUNCT__
#define __FUNCT__ "VertexCreate"
PetscErrorCode VertexCreate(FiberField field, Vertex *vertex)
{
  Vertex v;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = ArrayAppend( field->fibers, &v); CHKERRQ(ierr);
  ierr = UniqueIDGenerate(field->vid,&v->id); CHKERRQ(ierr);

  *vertex = v;
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "VertexDestroy"
PetscErrorCode VertexDestroy(FiberField field, Vertex v0 )
{
  int i;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  for ( i = 0; i < MAXEDGES; i++) {
    ierr = VertexRemoveEdge( v0, v0->v[i]); CHKERRQ(ierr);
  }

  // TODO: reset pointers when last vertex moved
  //ierr = ArrayDelete1( field->fibers, v0->idx ); CHKERRQ(ierr);
  
  SETERRQ(PETSC_COMM_SELF,0,"NOT IMP");

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

#undef __FUNCT__
#define __FUNCT__ "VertexRemoveEdge"
PetscErrorCode VertexRemoveEdge( Vertex v0, Vertex v1 )
{
  int j;
  //PetscErrorCode ierr;

  PetscFunctionBegin;
  if (v1 != NULL) {
    for (j = 0; j < MAXEDGES; j++) {
      v1->v[j] = NULL;
    }
  }
  PetscFunctionReturn(0);
}
