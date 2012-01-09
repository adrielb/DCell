#include "LevelSetMethod.h"
#include "LSM_private.h"

// Min-heap
static inline PetscReal FMMComparator( void *parent, void *child ) {
  return ((FMMNode)parent)->phi - ((FMMNode)child)->phi;
}

#undef __FUNCT__
#define __FUNCT__ "LevelSetReinitialize"
PetscErrorCode LevelSetReinitialize( LevelSet ls )
{
  static MemCache mc;
  static Heap heap;
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  ierr = PetscLogEventBegin(EVENT_LevelSetReinitialize,0,0,0,0); CHKERRQ(ierr);
  if( mc == NULL ) {
    ierr = MemCacheCreate(sizeof(struct _FMMNode),1e5,&mc); CHKERRQ(ierr);
    ierr = HeapCreate(FMMComparator, &heap); CHKERRQ(ierr);
  }
  
  // Ensure narrow band does not extend beyond allocated memory
  ierr = LevelSetResize(ls); CHKERRQ(ierr);

  int i;
  for( i = 0; i < ls->phi->SIZE; i++ ){
    ls->phi->v1[i] = ls->phi->v1[i] > 0 ? ls->PHI_INF : -ls->PHI_INF;
  }

  // The FMM will rebuild the narrow band
  ierr = ArraySetSize(ls->band,0); CHKERRQ(ierr);

  ierr = FMM_InitializeEikonal( ls, mc, heap ); CHKERRQ(ierr);
  
  // TODO: sort band according to x, y for better cache locality since FMM sorts band by phi
  ierr = PetscLogEventEnd(EVENT_LevelSetReinitialize,0,0,0,0); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "FMM_InitializeEikonal"
PetscErrorCode FMM_InitializeEikonal( LevelSet ls, MemCache mc, Heap heap )
{
  int i;
  int sign;
  double **phi2D, ***phi3D;
  int len = ArrayLength(ls->irregularNodes);
  IrregularNode *nodes = ArrayGetData(ls->irregularNodes);
  IrregularNode *n;
  iCoor *band;
  const PetscBool is2D = ls->phi->is2D;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = GridGet(ls->phi,&phi2D); CHKERRQ(ierr);
  ierr = GridGet(ls->phi,&phi3D); CHKERRQ(ierr);
//  Set the boundary condition for the fast marching method
  for( i = 0; i < len; i++ ) {
    n = &nodes[i];
    if( n->axis != -1 ) continue;
    ierr = ArrayAppend(ls->band,&band); CHKERRQ(ierr);
    band->x = n->pos.x;
    band->y = n->pos.y;
    band->z = n->pos.z;
    if( is2D ) {
      phi2D[n->pos.y][n->pos.x]          = n->signCenter * sqrt( PetscSqr(n->op.x) + PetscSqr(n->op.y) );
    } else {
      phi3D[n->pos.z][n->pos.y][n->pos.x]= n->signCenter * sqrt( PetscSqr(n->op.x) + PetscSqr(n->op.y) + PetscSqr(n->op.z));
    }
  }
// Add neighbors of BC's to heap
  for ( sign = -1; sign <= 1; sign+=2) {
    for( i = 0; i < len; i++ ) {
      n = &nodes[i];
       // if stencil node, skip it; (must be an ortho-proj node)
      if( n->axis != -1 ) continue;
      if( is2D ) {
        if( sign*phi2D[n->pos.y][n->pos.x] >= 0. ) {
          ierr = FMM_PushNeighbors2D(ls, mc, heap, sign, n->pos); CHKERRQ(ierr);
        }
      } else {
        if( sign*phi3D[n->pos.z][n->pos.y][n->pos.x] >= 0. ) {
          ierr = FMM_PushNeighbors3D(ls, mc, heap, sign, n->pos); CHKERRQ(ierr);
        }
      } // if bc same sign as domain
    } // for each irregular node i
    ierr = FMM_SolveEikonal( ls, mc, heap, sign ); CHKERRQ(ierr);
  } // for each sign=[-1,+1]

  PetscFunctionReturn(0);
}

// keep track of heap size for printing to PetscInfo
static int maxHeapSize;

#undef __FUNCT__
#define __FUNCT__ "FMM_SolveEikonal"
PetscErrorCode FMM_SolveEikonal( LevelSet ls, MemCache mc, Heap heap, int sign )
{
  PetscReal phi_val;
  PetscReal **phi2D, ***phi3D;
  FMMNode node;
  iCoor *band, pos;
  const PetscBool is2D = ls->phi->is2D;
  PetscErrorCode ierr=0;

  PetscFunctionBegin;
  maxHeapSize = 0;
  ierr = GridGet(ls->phi,&phi2D); CHKERRQ(ierr);
  ierr = GridGet(ls->phi,&phi3D); CHKERRQ(ierr);

  while( !HeapIsEmpty(heap) ) {
    ierr = HeapPop( heap, &node ); CHKERRQ(ierr);
    pos = node->pos;
    phi_val = node->phi;
    ierr = MemCacheFree(mc, node); CHKERRQ(ierr);
    if( phi_val > ls->bandWidth ) break;
    if( is2D ) {
      if( phi_val >= sign*phi2D[pos.y][pos.x] ) continue;
    } else {
      if( phi_val >= sign*phi3D[pos.z][pos.y][pos.x] ) continue;
    }
    // Add node to band (aka travel time is fixed)
    ierr = ArrayAppend(ls->band,&band); CHKERRQ(ierr);
    band->x = pos.x;
    band->y = pos.y;
    band->z = pos.z;

    if( is2D ) {
      phi2D[pos.y][pos.x] = sign * phi_val;
      ierr = FMM_PushNeighbors2D(ls, mc, heap, sign, pos); CHKERRQ(ierr);
    } else {
      phi3D[pos.z][pos.y][pos.x] = sign * phi_val;
      ierr = FMM_PushNeighbors3D(ls, mc, heap, sign, pos); CHKERRQ(ierr);
    }
  } // while within bandwidth

  ierr = PetscInfo1(0,"Max heap size = %d\n",maxHeapSize); CHKERRQ(ierr);

  // Dispose of all allocated node remaining on heap
  ierr = HeapClear(heap,mc); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}
