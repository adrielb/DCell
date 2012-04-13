#include "Common.h"

typedef struct _Quadtree *Quadtree;

struct _Quadtree {
  Quadtree parent;
  Quadtree quads[4];
  int level;
  Coor a,b; // [a.x---b.x]
  Coor p;   // [....p....]
  MemCache mem;
  void *data;
};

PetscErrorCode CreateParentQuad( Quadtree q, iCoor p, Quadtree *newquad );
PetscErrorCode CreateChildQuad(  Quadtree q, iCoor p, Quadtree *newquad );

inline PetscTruth Quadtree_InQuad( Quadtree q, iCoor p )
{
  return q->a.x <= p.x && p.x < q->b.x &&
         q->a.y <= p.y && p.y < q->b.y;
}

#undef __FUNCT__
#define __FUNCT__ "QuadtreeCreate"
PetscErrorCode QuadtreeCreate( int maxDepth, iCoor a, iCoor b, Quadtree *root )
{
  PetscErrorCode ierr;

  PetscFunctionBegin;

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "QuadtreeDestroy"
PetscErrorCode QuadtreeDestroy( int maxDepth, iCoor a, iCoor b, Quadtree *root )
{
  PetscErrorCode ierr;

  PetscFunctionBegin;

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "QuadtreeInsert"
PetscErrorCode QuadtreeInsert( Quadtree q, iCoor p, void **data )
{
  int i;
  PetscErrorCode ierr;

  PetscFunctionBegin;

  do {
    if( Quadtree_InQuad( q, p ) ) {
      if( q->level == 0 ) {
        *data = q->data;
        break;
      }
      // p is in quad, but not at level 0
      for (i = 0; i < 4; ++i) {
        if( q->quads[i] == NULL ) {
          ierr = CreateChildQuad(q,p,&q->quads[i]); CHKERRQ(ierr);
          q = q->quads[i];
          break;
        }
        if( Quadtree_InQuad( q->quads[i], p) ) {
          q = q->quads[i];
          break;
        }
      }
    } else { // not in current quad
      // go up to parent
      if( q->parent == NULL ) {
        ierr = CreateParentQuad(q,p,&q->parent); CHKERRQ(ierr);
      }
      q = q->parent;
    }
  } while (PETSC_TRUE);

  PetscFunctionReturn(0);
}

PetscErrorCode CreateParentQuad(Quadtree q, iCoor p, Quadtree *newquad ) {
  Quadtree qp;
  PetscErrorCode ierr;
  ierr = MemCacheAlloc(q->mem, &qp); CHKERRQ(ierr);
  qp->mem = q->mem;
  qp->quads[0] = q;
  qp->level = q->level + 1;
  if( p.x < q->p.x ) { // [-------*------ax-----px-----bx]
    qp->a.x = 2 * q->a.x - q->b.x;
    qp->b.x = q->b.x;
  } else {             // [ax-----px-----bx-------*------]
    qp->a.x = q->a.x;
    qp->b.x = 2 * q->b.x - q->a.x;
  }
  if( p.y < q->p.y ) { // [-------*------ax-----px-----bx]
    qp->a.y = 2 * q->a.y - q->b.y;
    qp->b.x = q->b.y;
  } else {             // [ax-----px-----bx-------*------]
    qp->a.y = q->a.y;
    qp->b.y = 2 * q->b.y - q->a.y;
  }
  qp->p.x = ( qp->a.x + qp->b.x ) / 2;
  qp->p.y = ( qp->a.y + qp->b.y ) / 2;
  *newquad = qp;
}

PetscErrorCode CreateChildQuad(Quadtree q, iCoor p, Quadtree *newquad ) {
  Quadtree qi;
  PetscErrorCode ierr;
  ierr = MemCacheAlloc(q->mem,(void*)&qi); CHKERRQ(ierr);
  qi->mem = q->mem;
  qi->parent = q;
  qi->level = q->level - 1;
  if( p.x < q->p.x ) { // [ax--*--px]-----bx
    qi->a.x = q->a.x;
    qi->b.x = q->p.x;
  } else {             // ax-----[px--*--bx]
    qi->a.x = q->p.x;
    qi->b.x = q->b.x;
  }
  if( p.y < q->p.y ) {
    qi->a.y = q->a.y;
    qi->b.y = q->p.y;
  } else {
    qi->a.y = q->p.y;
    qi->b.y = q->b.y;
  }
  qi->p.x = ( qi->a.x + qi->b.x ) / 2;
  qi->p.y = ( qi->a.y + qi->b.y ) / 2;
}

#undef __FUNCT__
#define __FUNCT__ "QuadtreeGet"
PetscErrorCode QuadtreeGet( Quadtree q, iCoor p, void **data )
{
  int i;
  do {
    if( Quadtree_InQuad( q, p ) ) {
      if( q->level == 0 ) {
        *data = q->data;
        break;
      }
      for (i = 0; i < 4; ++i) {
        if( Quadtree_InQuad( q->quads[i], p) ) {
          q = q->quads[i];
          break;
        }
      }

    } else { // not in current quad
      // go up to parent
      q = q->parent;
      if( q == NULL ) {
        // Error: requested elem outside root of tree
        break;
      }
    }
  } while (PETSC_TRUE);
}


#undef __FUNCT__
#define __FUNCT__ "QuadtreeGetSphere"
PetscErrorCode QuadtreeGetSphere( Quadtree root, Coor center, PetscReal radius, void **elems, int *len )
{
  PetscErrorCode ierr;

  PetscFunctionBegin;

  PetscFunctionReturn(0);
}
