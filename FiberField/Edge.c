#include "FiberField.h"

#undef __FUNCT__
#define __FUNCT__ "EdgeCreate"
PetscErrorCode EdgeCreate( FiberField fibers, Edge *edge  )
{
  PetscErrorCode ierr;

  PetscFunctionBegin;

  ierr = MemCacheAlloc(fibers->mcEdges,&e); CHKERRQ(ierr);
  ierr = UniqueIDGenerate(fibers->eid,&e->id); CHKERRQ(ierr);
  e->a = v1;
  e->id_a = v1->id;
  e->b = v2;
  e->id_b = v2->id;
  ierr = VertexAddEdge(v1,e); CHKERRQ(ierr);
  ierr = VertexAddEdge(v2,e); CHKERRQ(ierr);

  *edge = e;
  PetscFunctionReturn(0);
}
