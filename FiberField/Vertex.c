#include "FiberField.h"

PetscErrorCode Vertex_Link( Vertex v0, Edge e ); 

#undef __FUNCT__
#define __FUNCT__ "FiberFieldAddVertex"
PetscErrorCode FiberFieldAddVertex(FiberField field, FiberTypeID vtype, Vertex *vertex)
{
  int i;
  Vertex v;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = ArrayAppend( field->verts, &v); CHKERRQ(ierr);
  ierr = UniqueIDGenerate(field->vid,&v->vID); CHKERRQ(ierr);

  for (i = 0; i < MAXEDGES; i++) {
    v->eID[i] = FIBERFIELD_NO_EDGE;
  }

  FiberType *ftype;
  ierr = ArrayGet( field->fibertypesDB, vtype, &ftype); CHKERRQ(ierr);
  if (ftype->isEdge) {
    SETERRQ( PETSC_COMM_SELF, 0, "Wrong fiber type: adding vertex with ftype edge\n");
  }

  v->type = vtype;

  *vertex = v;
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "FiberFieldAddEdge"
PetscErrorCode FiberFieldAddEdge( FiberField f,  Vertex v0, Vertex v1, FiberTypeID etype, PetscReal l0 )
{
  Edge edge;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  FiberType *ftype;
  ierr = ArrayGet(f->fibertypesDB, etype, &ftype); CHKERRQ(ierr);
  if (!ftype->isEdge) {
    SETERRQ( PETSC_COMM_SELF, 0, "Wrong fiber type: adding edge with ftype vertex\n");
  }

  ierr = ArrayAppend( f->edges, &edge); CHKERRQ(ierr);
  ierr = UniqueIDGenerate( f->eid, &edge->eID); CHKERRQ(ierr); 
  edge->vID[0] = v0->vID;
  edge->vID[1] = v1->vID;
  edge->type = etype;
  edge->l0 = l0;

  ierr = Vertex_Link( v0, edge ); CHKERRQ(ierr);
  ierr = Vertex_Link( v1, edge ); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "Vertex_Link"
PetscErrorCode Vertex_Link( Vertex v0, Edge e )
{
  int i;
  //PetscErrorCode ierr;

  PetscFunctionBegin;

  // Add edge into v0's edge list
  for (i = 0; i < MAXEDGES; ++i) {
    if( v0->eID[i] == FIBERFIELD_NO_EDGE ) {
      v0->eID[i] = e->eID;
      break;
    }
  }
  if( i == MAXEDGES ) {
    SETERRQ1(PETSC_COMM_SELF, 0, "MAXEDGES (%d) reached\n", MAXEDGES);
  }
  PetscFunctionReturn(0);
}

#define __BUGGY__
#ifndef __BUGGY__ 


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

#endif
