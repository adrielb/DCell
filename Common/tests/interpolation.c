#include "Common.h"
#include "Grid.h"

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc, char **args)
{
  PetscErrorCode  ierr;
  PetscFunctionBegin;
  ierr = DCellInit(); CHKERRQ(ierr);

  Grid g;
  Coor dh = {1.3, 1.3, 0.0};
  iCoor pos  = {0, 0, 0};
  iCoor size = {4, 4, 0};
  GridCreate( dh, pos, size, 1, &g);
  PetscReal **h;
  GridGet(g, &h);
  h[0][0] =-1; h[0][1] = -1; h[0][2] =-1; h[0][3] =-1;
  h[1][0] =-1; h[1][1] =  1; h[1][2] = 1; h[1][3] =-1;
  h[2][0] =-1; h[2][1] = -1; h[2][2] = 2; h[2][3] =-1;
  h[3][0] =-1; h[3][1] = -1; h[3][2] =-1; h[3][3] =-1;

  PetscReal dI = 0.1;
  Coor di = {dI, dI, dI};
  AABB aabb = g->aabb;
  aabb.lo.x -= 1;
  aabb.lo.y -= 1;
  aabb.hi.x += 1;
  aabb.hi.y += 1;
  CoorToIndex(aabb.lo, di, aabb.hi, &size);
  Grid gi;
  GridCreate( di, (iCoor){0,0,0}, size, 1, &gi);
  PetscReal **a;
  GridGet(gi, &a);

  g->Interpolate = GridInterpolate_Cubic;

  Coor X;
  iCoor x,p,q;
  GridGetBounds(gi, &p, &q);
  for (x.y = p.y; x.y < q.y; ++x.y) {
    for (x.x = p.x; x.x < q.x; ++x.x) {
      IndexToCoor(aabb.lo,di,x,&X);
      ierr = GridInterpolate( g, X, &a[x.y][x.x]); CHKERRQ(ierr);
    }
  }
  GridWrite( gi, 0);


  ierr = GridDestroy(gi); CHKERRQ(ierr);
  ierr = GridDestroy(g); CHKERRQ(ierr);
  ierr = DCellFinalize(); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
