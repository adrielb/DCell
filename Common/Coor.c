#include "Common.h"

inline void CoorToIndex( const Coor origin, const Coor dh, const Coor point, iCoor* P )
{
  P->x = (int)floor( (point.x - origin.x) / dh.x );
  P->y = (int)floor( (point.y - origin.y) / dh.y );
  P->z = (int)floor( (point.z - origin.z) / dh.z );
}

inline void CoorToIndex2( const Coor origin, const Coor dh, const Coor point, Coor* p, iCoor* P)
{
  p->x = (point.x - origin.x) / dh.x;
  p->y = (point.y - origin.y) / dh.y;
  p->z = (point.z - origin.z) / dh.z;

  P->x = (int)floor( p->x );
  P->y = (int)floor( p->y );
  P->z = (int)floor( p->z );
}

inline PetscBool AABBPointInBox( const AABB box, const Coor p )
{
  return  p.x < box.lo.x || box.hi.x < p.x ||
          p.y < box.lo.y || box.hi.y < p.y ||
          p.z < box.lo.z || box.hi.z < p.z;
}
