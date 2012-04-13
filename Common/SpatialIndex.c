#include "Common.h"

struct _SpatialIndex {
  Coor lo, hi; // [lo.x---hi.x]
  Coor dh;     // [lo.x---lo.x+dh.x---lo.x+2*dh.x...]
  iCoor numbins;  // [3x4x5]
  iCoor a,b;   // start bin coor, end bin coor
  Array bins;  // 2D/3D array of bins
  MemCache items;
};

struct _SpatialItem {
  SpatialItem next;
  AABB box;
  Coor p;
  void *item;
};

PetscErrorCode SpatialIndexBin( SpatialIndex sidx, Coor pt, iCoor *bin );

#undef __FUNCT__
#define __FUNCT__ "SpatialIndexCreate"
PetscErrorCode SpatialIndexCreate( const char name[], SpatialIndex *sidx )
{
  size_t chunksize = 65536; // TODO: make petsc option for this
  SpatialIndex s;
  char tmp[256];
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscNew(struct _SpatialIndex, &s); CHKERRQ(ierr);
  sprintf(tmp, "%s_%s", name, "spatialbins");
  ierr = ArrayCreate( tmp, sizeof(SpatialItem), &s->bins); CHKERRQ(ierr);
  sprintf(tmp, "%s_%s", name, "spatialItems");
  ierr = MemCacheCreate( tmp, sizeof(struct _SpatialItem), chunksize, &s->items); CHKERRQ(ierr);
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
  s->numbins.x = ( hi.x - lo.x ) / dh.x + 2;
  s->numbins.y = ( hi.y - lo.y ) / dh.y + 2;
  s->numbins.z = ( hi.z - lo.z ) / dh.z + 2;
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
  ierr = ArrayDestroy( sidx->bins ); CHKERRQ(ierr);
  ierr = MemCacheDestroy( sidx->items ); CHKERRQ(ierr);
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
PetscErrorCode SpatialIndexInsertPoint( SpatialIndex sidx, Coor pt, void *data )
{
  iCoor bin;
  SpatialItem *iter, newitem;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = SpatialIndexBin( sidx, pt, &bin ); CHKERRQ(ierr);

  ierr = ArrayGetCoor( sidx->bins, bin, &iter); CHKERRQ(ierr);

  ierr = MemCacheAlloc(sidx->items,&newitem); CHKERRQ(ierr);
  newitem->next = PETSC_NULL;
  newitem->p = pt;
  newitem->item = data;

  // Prepend to linked list
  if( *iter != PETSC_NULL ) {
    newitem->next = *iter;
  }

  *iter = newitem;

  PetscFunctionReturn(0);
}

PetscErrorCode SpatialIndexBin( SpatialIndex sidx, Coor pt, iCoor *bin )
{
  bin->x = ( pt.x - sidx->lo.x ) / sidx->dh.x;
  bin->y = ( pt.y - sidx->lo.y ) / sidx->dh.y;
  bin->z = ( pt.z - sidx->lo.z ) / sidx->dh.z;

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
PetscErrorCode SpatialIndexQueryPoints( SpatialIndex sidx, Coor center, PetscReal radius, const int MAXLEN, int *len, void *items[] )
{
  int c = 0;
  iCoor lo, hi, b;
  Coor p;
  PetscReal r2 = radius * radius;
  PetscReal sqdist;
  SpatialItem iter;

  PetscErrorCode ierr;

  PetscFunctionBegin;
  lo.x = ( center.x-radius - sidx->lo.x ) / sidx->dh.x;
  lo.y = ( center.y-radius - sidx->lo.y ) / sidx->dh.y;
  lo.z = ( center.z-radius - sidx->lo.z ) / sidx->dh.z;
  if( lo.x < sidx->a.x ) lo.x = sidx->a.x;
  if( lo.y < sidx->a.y ) lo.y = sidx->a.y;
  if( lo.z < sidx->a.z ) lo.z = sidx->a.z;

  hi.x = ( center.x+radius - sidx->lo.x ) / sidx->dh.x;
  hi.y = ( center.y+radius - sidx->lo.y ) / sidx->dh.y;
  hi.z = ( center.z+radius - sidx->lo.z ) / sidx->dh.z;
  if( sidx->b.x < hi.x ) hi.x = sidx->b.x;
  if( sidx->b.y < hi.y ) hi.y = sidx->b.y;
  if( sidx->b.z < hi.z ) hi.z = sidx->b.z;

  for (b.z = lo.z; b.z <= hi.z; ++b.z) {
    for (b.y = lo.y; b.y <= hi.y; ++b.y) {
      for (b.x = lo.x; b.x <= hi.x; ++b.x) {
        ierr = ArrayGetCoorP( sidx->bins, b, &iter); CHKERRQ(ierr);
        while( iter != PETSC_NULL ) {
          p = iter->p;
          p.x = p.x - center.x;
          p.y = p.y - center.y;
          p.z = p.z - center.z;
          sqdist = p.x*p.x + p.y*p.y + p.z*p.z;
          if( sqdist <= r2 ) {
            items[c] = iter->item;
            c++;
            if( c == MAXLEN ) {
              ierr = PetscInfo1(0,"Max query size reached: MAXLEN = %d\n", MAXLEN); CHKERRQ(ierr);
              *len = c;
              PetscFunctionReturn(0);
            } // c == Np
          } // sqdist <= r2
          iter = iter->next;
        } // iter
      } // x
    } // y
  } // z

  *len = c;

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "SpatialIndexClear"
PetscErrorCode SpatialIndexClear( SpatialIndex sidx )
{
  int i;
  int len = ArrayLength(sidx->bins);
  SpatialItem item, nextitem, *items = ArrayGetData(sidx->bins);
  PetscErrorCode ierr;

  PetscFunctionBegin;
  for (i = 0; i < len; ++i) {
    if( items[i] != PETSC_NULL ) {
      item = items[i];
      while( item != PETSC_NULL ) {
        nextitem = item->next;
        ierr = MemCacheFree( sidx->items, item); CHKERRQ(ierr);
        item = nextitem;
      }
    }
  }
  ierr = ArrayZero(sidx->bins); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "SpatialIndexCollide"
PetscErrorCode SpatialIndexCollide( SpatialIndex sidx, AABB box, SpatialItem items )
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
