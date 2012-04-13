#include "LevelSetMethod.h"

#undef __FUNCT__
#define __FUNCT__ "LevelSetResize"
PetscErrorCode LevelSetResize(LevelSet ls)
{
  int i,j,k;
  iCoor max, min;
  IrregularNode *n;
  iCoor shift;    // shift amount in new memory
  iCoor dim;      // minimum size of new memory to hold ls
  iCoor p,q;  // Current size of grid
  const int thres = ls->bandWidth + 1; // threshold that bounding box is on border
  const int BUF = thres + 5; // 5 grids of additional buffer on each edge
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = GridGetBounds(ls->phi,&p,&q); CHKERRQ(ierr);

  // Find bounding box {min, max} of irregular nodes
  max = p;
  min = q;
  for ( i = 0; i < ArrayLength(ls->irregularNodes); ++i) {
    ierr = ArrayGet(ls->irregularNodes,i,(void*)&n); CHKERRQ(ierr);
    if( max.x < n->pos.x ) max.x = n->pos.x;
    if( min.x > n->pos.x ) min.x = n->pos.x;
    if( max.y < n->pos.y ) max.y = n->pos.y;
    if( min.y > n->pos.y ) min.y = n->pos.y;
    if( max.z < n->pos.z ) max.z = n->pos.z;
    if( min.z > n->pos.z ) min.z = n->pos.z;
  }

  // If bb does not touch boundary, quit, nothing to do
  if( p.x <= min.x - thres && max.x + thres < q.x &&
      p.y <= min.y - thres && max.y + thres < q.y ) {
    if( ls->phi->is2D ) {
      PetscFunctionReturn(0);
    } else {
      if( p.z <= min.z - thres && max.z + thres < q.z ) {
        PetscFunctionReturn(0);
      }
    }
  }

  /* else
   * bounding box touches boundary
   * 1) Resize the memory extents
   * 2) Shift LS to center of new memory
   */

  dim.x = max.x - min.x + 1 + 2*BUF;
  dim.y = max.y - min.y + 1 + 2*BUF;
  dim.z = max.z - min.z + 1 + 2*BUF;

  shift.x = min.x - BUF;
  shift.y = min.y - BUF;
  shift.z = min.z - BUF;

  // Shift and resize tmp grid
  ierr = GridResize(ls->tmp,shift,dim); CHKERRQ(ierr);

  // Fill +PHI_INF for grid points outside narrow band
  ierr = VecSet(ls->tmp->v,ls->PHI_INF); CHKERRQ(ierr);

  // Calculate the union of the old and new extents
  p.x = PetscMax(p.x,shift.x);
  p.y = PetscMax(p.y,shift.y);
  p.z = PetscMax(p.z,shift.z);
  q.x = PetscMin(q.x,shift.x+dim.x);
  q.y = PetscMin(q.y,shift.y+dim.y);
  q.z = PetscMin(q.z,shift.z+dim.z);

  // Copy old ls into new mem tmp
  if( ls->phi->is2D ) {
    PetscReal **phi, **tmp;
    ierr = GridGet(ls->phi,&phi); CHKERRQ(ierr);
    ierr = GridGet(ls->tmp,&tmp); CHKERRQ(ierr);
    for (j = p.y; j < q.y; ++j) {
      for (i = p.x; i < q.x; ++i) {
        tmp[j][i] = phi[j][i];
      }
    }
    dim.z = 0;
  } else {
    PetscReal ***phi, ***tmp;
    ierr = GridGet(ls->phi,&phi); CHKERRQ(ierr);
    ierr = GridGet(ls->tmp,&tmp); CHKERRQ(ierr);
    for (k = p.z; k < q.z; ++k) {
      for (j = p.y; j < q.y; ++j) {
        for (i = p.x; i < q.x; ++i) {
          tmp[k][j][i] = phi[k][j][i];
        }
      }
    }
  }

  // Copy level set from tmp to phi
  ierr = GridResize(ls->phi,shift,dim); CHKERRQ(ierr);
  ierr = VecCopy(ls->tmp->v,ls->phi->v); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
