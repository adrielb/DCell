#include "FiberField.h"

#undef __FUNCT__
#define __FUNCT__ "VertexCreate"
PetscErrorCode VertexCreate(FiberField field, Vertex *vertex)
{
  Vertex v;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = ArrayAppend( field->fibers, &v); CHKERRQ(ierr);
  ierr = UniqueIDGenerate(field->vid,&v->ID); CHKERRQ(ierr);

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
  PetscErrorCode ierr;
  PetscFunctionBegin;
  ierr = Vertex_Link( v0, v1, etype); CHKERRQ(ierr);
  ierr = Vertex_Link( v1, v0, etype); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "Vertex_Link"
PetscErrorCode Vertex_Link( Vertex v0, Vertex v1, EdgeType etype )
{
  int i;
  //PetscErrorCode ierr;

  PetscFunctionBegin;

  // Add v1 into v0's vertex list
  for (i = 0; i < MAXEDGES; ++i) {
    if( v0->v[i] == NULL ) {
      v0->v[i] = v1;
      v0->vID[i] = v1->ID;
      v0->e[i] = etype;
      break;
    }
  }
  if( i == MAXEDGES ) {
    SETERRQ1(PETSC_COMM_SELF, 0, "MAXEDGES (%d) reached\n", MAXEDGES);
  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "VertexRemoveEdge"
PetscErrorCode VertexRemoveEdge( Vertex v0, Vertex v1 )
{
  PetscErrorCode ierr;

  PetscFunctionBegin;

  ierr = Vertex_Unlink( v0, v1->ID); CHKERRQ(ierr);
  ierr = Vertex_Unlink( v1, v0->ID); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "Vertex_Unlink"
PetscErrorCode Vertex_Unlink( Vertex v0, VertexID v1_ID )
{
  int j;
  //PetscErrorCode ierr;

  PetscFunctionBegin;

  // remove v1 from v0's vertex list
  for (j = 0; j < MAXEDGES; j++) {
    if (v0->vID[j] == v1_ID) {
      v0->v[j] = NULL;
      v0->vID[j] = FIBERFIELD_NO_VERTEX;
      v0->e[j] = FIBERFIELD_NO_VERTEX;
      break;
    }
  }
  if (j == MAXEDGES) {
    SETERRQ2( PETSC_COMM_SELF, 0, "Verticies not connected: v0 = %d, v1 = %d\n", v0->ID, v1_ID );
  }
  
  PetscFunctionReturn(0);
}

