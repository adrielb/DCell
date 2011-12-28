#include "Common.h"

typedef struct _SpatialIndex *SpatialIndex;

struct _SpatialIndex {
  Coor lo, hi; // [lo.x---hi.x]
  Coor dh;     // [lo.x---lo.x+dh.x---lo.x+2*dh.x...]
  void *grid;  // 2D/3D array of bins

  MemCache mem;
};

typedef struct _SpatialItem {
  SpatialItem next;
  AABB box;
  void *item
} *SpatialItem;


#undef __FUNCT__
#define __FUNCT__ "SpatialIndexCreate"
PetscErrorCode SpatialIndexCreate( Coor lo, Coor hi, Coor dh, SpatialIndex *sidx )
{
  PetscErrorCode ierr;

  PetscFunctionBegin;


  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "SpatialIndexDestroy"
PetscErrorCode SpatialIndexDestroy( SpatialIndex sidx )
{
  PetscErrorCode ierr;

  PetscFunctionBegin;


  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "SpatialIndexInsertBox"
PetscErrorCode SpatialIndexInsertBox( SpatialIndex sidx, AABB box, void *item )
{
  Coor center;
  iCoor bin;
  PetscErrorCode ierr;

  PetscFunctionBegin;

  center.x = box.hi.x - box.lo.x;

  bin.x = ( center.x - sidx->lo.x ) / sidx->dh.x;

  newitem = grid[bin.x];
  grid[bin.x] = newitem;;

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "SpatialIndexCollide"
PetscErrorCode SpatialIndexCollide( SpatialIndex sidx, AABB box, SpatialBoxItem items )
{
  PetscErrorCode ierr;

  PetscFunctionBegin;


  PetscFunctionReturn(0);
}
