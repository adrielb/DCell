#include "Common.h"

inline void CoorToIndex( const Coor origin, const Coor dh, const Coor point, iCoor *P )
{
  P->x = (int)floor( (point.x - origin.x) / dh.x );
  P->y = (int)floor( (point.y - origin.y) / dh.y );
  P->z = (int)floor( (point.z - origin.z) / dh.z );
}

inline void CoorToIndex2( const Coor origin, const Coor dh, const Coor point, iCoor* P, Coor *p )
{
  p->x = (point.x - origin.x) / dh.x;
  p->y = (point.y - origin.y) / dh.y;
  p->z = (point.z - origin.z) / dh.z;

  P->x = (int)floor( p->x );
  P->y = (int)floor( p->y );
  P->z = (int)floor( p->z );
}

inline void IndexToCoor( const Coor origin, const Coor dh, const iCoor i, Coor *p )
{
  p->x = origin.x + i.x * dh.x;
  p->y = origin.y + i.y * dh.y;
  p->z = origin.z + i.z * dh.z;
}

inline PetscBool AABBPointInBox( const AABB box, const Coor p )
{
  const PetscBool notInBox =
      p.x < box.lo.x || box.hi.x < p.x ||
      p.y < box.lo.y || box.hi.y < p.y ||
      p.z < box.lo.z || box.hi.z < p.z;
  return !notInBox;
}

inline PetscReal CoorGet( const Coor x, const int i )
{
  return ((PetscReal*)&x)[i];
}

inline void CoorSet( const Coor *x, const int i, const PetscReal val )
{
  ((PetscReal*)x)[i] = val;
}

inline PetscInt iCoorGet( const iCoor x, const int i )
{
  return ((PetscInt*)&x)[i];
}
