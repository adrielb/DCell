#include "Common.h"

struct _SpatialIndex {
  Coor lo, hi; // [lo.x---hi.x]
  Coor dh;     // [lo.x---lo.x+dh.x---lo.x+2*dh.x...]
  iCoor numbins;  // [3x4x5]
  iCoor a,b;   // start bin coor, end bin coor
  Array bins;  // 2D/3D array of bins
  Array queriedItems;
  size_t offset;
};

struct _SpatialItem {
  SpatialItem next;
  AABB *box;
  Coor *p;
};

PetscErrorCode SpatialIndexBin( SpatialIndex sidx, Coor pt, iCoor *bin );
PetscErrorCode SpatialIndex_AddItem( SpatialIndex sidx, void *newItem );

#undef __FUNCT__
#define __FUNCT__ "SpatialIndexCreate"
PetscErrorCode SpatialIndexCreate( const char name[], SpatialIndex *sidx )
{
  SpatialIndex s;
  char tmp[256];
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscNew(struct _SpatialIndex, &s); CHKERRQ(ierr);
  sprintf(tmp, "%s_%s", name, "spatialbins");
  ierr = ArrayCreate( tmp, sizeof(SpatialItem), &s->bins); CHKERRQ(ierr);
  sprintf(tmp, "%s_%s", name, "spatialItems");
  ierr = ArrayCreate( "queriedItems", sizeof(void*), &s->queriedItems ); CHKERRQ(ierr);
  ierr = SpatialIndexSetDomain( s, (Coor){0,0,0}, (Coor){1,1,1}, (Coor){1,1,1}); CHKERRQ(ierr);
  *sidx = s;
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "SpatialIndexSetDomain"
PetscErrorCode SpatialIndexSetDomain( SpatialIndex s, Coor lo, Coor hi, Coor dh )
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  s->lo = lo;
  s->hi = hi;
  s->dh = dh;
  CoorToIndex( lo, dh, hi, &s->numbins );
  s->numbins.x += 2;
  s->numbins.y += 2;
  s->numbins.z += 2;
  s->a = (iCoor){-1,-1,-1};
  s->b.x = s->numbins.x - 2;
  s->b.y = s->numbins.y - 2;
  s->b.z = s->numbins.z - 2;
  ierr = ArraySetCoor( s->bins, s->a, s->numbins); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "SpatialIndexDestroy"
PetscErrorCode SpatialIndexDestroy( SpatialIndex sidx )
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = ArrayDestroy( sidx->queriedItems ); CHKERRQ(ierr);
  ierr = ArrayDestroy( sidx->bins ); CHKERRQ(ierr);
  ierr = PetscFree( sidx ); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "SpatialIndexInsertBox"
PetscErrorCode SpatialIndexInsertBox( SpatialIndex sidx, AABB box, void *item )
{
/*
  Coor center;
  iCoor bin;
  PetscErrorCode ierr;

  PetscFunctionBegin;

  center.x = box.hi.x - box.lo.x;

  ierr = SpatialIndexBin( sidx, center, &bin ); CHKERRQ(ierr);

  newitem = grid[bin.x];
  grid[bin.x] = newitem;;
*/
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "SpatialIndexInsertPoint"
PetscErrorCode SpatialIndexInsertPoint( SpatialIndex sidx, Coor *pt, SpatialItem newitem )
{
  iCoor bin;
  SpatialItem *iter;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = SpatialIndexBin( sidx, *pt, &bin ); CHKERRQ(ierr);

  ierr = ArrayGetCoor( sidx->bins, bin, &iter); CHKERRQ(ierr);

  newitem->next = PETSC_NULL;
  newitem->p = pt;

  // Prepend to linked list
  if( *iter != PETSC_NULL ) {
    newitem->next = *iter;
  }

  *iter = newitem;

  PetscFunctionReturn(0);
}

PetscErrorCode SpatialIndexBin( SpatialIndex sidx, Coor pt, iCoor *bin )
{
  CoorToIndex( sidx->lo, sidx->dh, pt, bin );

  // if point outside domain, put in boundary bins
  if( bin->x < sidx->a.x ) bin->x = sidx->a.x;
  if( bin->y < sidx->a.y ) bin->y = sidx->a.y;
  if( bin->z < sidx->a.z ) bin->z = sidx->a.z;

  if( bin->x > sidx->b.x ) bin->x = sidx->b.x;
  if( bin->y > sidx->b.y ) bin->y = sidx->b.y;
  if( bin->z > sidx->b.z ) bin->z = sidx->b.z;

  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "SpatialIndexQueryPoints"
PetscErrorCode SpatialIndexQueryPoints( SpatialIndex sidx, Coor center, PetscReal radius, Array *items )
{
  iCoor lo, hi, b;
  Coor p;
  const PetscReal r2 = radius * radius;
  PetscReal sqdist;
  SpatialItem iter;
  PetscErrorCode ierr;

  PetscFunctionBegin;

//  CoorToIndex( sidx->lo.x, sidx->dh, center-radius, &lo );
  lo.x = (PetscInt)(( center.x-radius - sidx->lo.x ) / sidx->dh.x);
  lo.y = (PetscInt)(( center.y-radius - sidx->lo.y ) / sidx->dh.y);
  lo.z = (PetscInt)(( center.z-radius - sidx->lo.z ) / sidx->dh.z);
  if( lo.x < sidx->a.x ) lo.x = sidx->a.x;
  if( lo.y < sidx->a.y ) lo.y = sidx->a.y;
  if( lo.z < sidx->a.z ) lo.z = sidx->a.z;

  hi.x = (PetscInt)(( center.x+radius - sidx->lo.x ) / sidx->dh.x);
  hi.y = (PetscInt)(( center.y+radius - sidx->lo.y ) / sidx->dh.y);
  hi.z = (PetscInt)(( center.z+radius - sidx->lo.z ) / sidx->dh.z);
  if( sidx->b.x < hi.x ) hi.x = sidx->b.x;
  if( sidx->b.y < hi.y ) hi.y = sidx->b.y;
  if( sidx->b.z < hi.z ) hi.z = sidx->b.z;

  ierr = ArraySetSize(sidx->queriedItems, 0); CHKERRQ(ierr);
  for (b.z = lo.z; b.z <= hi.z; ++b.z) {
    for (b.y = lo.y; b.y <= hi.y; ++b.y) {
      for (b.x = lo.x; b.x <= hi.x; ++b.x) {
        ierr = ArrayGetCoorP( sidx->bins, b, &iter); CHKERRQ(ierr);
        while( iter != PETSC_NULL ) {
          p = *iter->p;
          p.x = p.x - center.x;
          p.y = p.y - center.y;
          p.z = p.z - center.z;
          sqdist = p.x*p.x + p.y*p.y + p.z*p.z;
          if( sqdist <= r2 ) {
            ierr = SpatialIndex_AddItem( sidx, iter ); CHKERRQ(ierr);
          } // sqdist <= r2
          iter = iter->next;
        } // iter
      } // x
    } // y
  } // z

  *items = sidx->queriedItems;

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "SpatialIndexQueryPointsBox"
PetscErrorCode SpatialIndexQueryPointsBox( SpatialIndex sidx, AABB box, Array *items )
{
  iCoor lo, hi, b;
  SpatialItem iter;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = SpatialIndexBin( sidx, box.lo, &lo ); CHKERRQ(ierr);
  ierr = SpatialIndexBin( sidx, box.hi, &hi ); CHKERRQ(ierr);

  ierr = ArraySetSize(sidx->queriedItems, 0); CHKERRQ(ierr);
  for (b.z = lo.z; b.z <= hi.z; ++b.z) {
    for (b.y = lo.y; b.y <= hi.y; ++b.y) {
      for (b.x = lo.x; b.x <= hi.x; ++b.x) {
        ierr = ArrayGetCoorP( sidx->bins, b, &iter); CHKERRQ(ierr);
        while( iter != PETSC_NULL ) {
          if( AABBPointInBox(box, *iter->p) ) {
            ierr = SpatialIndex_AddItem( sidx, iter ); CHKERRQ(ierr);
          }
          iter = iter->next;
        } // iter
      } // x
    } // y
  } // z

  *items = sidx->queriedItems;

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "SpatialIndexClear"
PetscErrorCode SpatialIndexClear( SpatialIndex sidx )
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = ArrayZero(sidx->bins); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "SpatialIndexCollide"
PetscErrorCode SpatialIndexCollide( SpatialIndex sidx, AABB box, void *items )
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = 0;
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "SpatialIndexPrint"
PetscErrorCode SpatialIndexPrint( SpatialIndex sidx )
{
  int len = 0;
  iCoor a = sidx->a;
  iCoor b = sidx->b;
  iCoor x;
  SpatialItem iter;

  int max = 0;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  for (x.z = a.z; x.z <= b.z; x.z+=1) {
    for (x.y = a.y; x.y <= b.y; x.y+=1) {
      for (x.x = a.x; x.x <= b.x; x.x+=1) {
        ierr = ArrayGetCoorP( sidx->bins, x, &iter); CHKERRQ(ierr);
        len = 0;
        while( iter != PETSC_NULL ) {
          len++;
          iter = iter->next;
        }
        if( len > max )
          max = len;
//          printf("%d %d %d %d\n", x.x, x.y, x.z, len);
      }
    }
  }
  printf("%d \n", max);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "SpatialIndex_AddItem"
PetscErrorCode SpatialIndex_AddItem( SpatialIndex sidx, void *newItem )
{
  void **item;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = ArrayAppend(sidx->queriedItems, &item); CHKERRQ(ierr);
  *item = (void*)((char*)newItem - sidx->offset);
  PetscFunctionReturn(0);
}
