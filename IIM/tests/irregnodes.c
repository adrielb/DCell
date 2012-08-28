#include "ImmersedInterfaceMethod.h"
#include "IIM_private.h"

PetscErrorCode findroots( Grid g, const Coor p, const Coor dr, Array roots );

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc, char **args)
{
  PetscErrorCode  ierr;
  PetscFunctionBegin;
  ierr = DCellInit(); CHKERRQ(ierr);

  const Coor dg = {1,1,1};
  const iCoor size = {3,3,0};
  const iCoor pos = {0,0,0};
  Grid g;
  ierr = GridCreate( dg, pos, size, 1, &g); CHKERRQ(ierr);
  PetscReal **val;
  ierr = GridGet(g, &val); CHKERRQ(ierr);
  val[2][0] = -1; val[2][1] = -1; val[2][2] =  1;
  val[1][0] = -1; val[1][1] = -1; val[1][2] =  1;
  val[0][0] = -1; val[0][1] =  1; val[0][2] =  1;
  ierr = GridWrite(g, 0); CHKERRQ(ierr);

  Array roots;
  ierr = ArrayCreate("roots", sizeof(Coor), &roots); CHKERRQ(ierr);

  const Coor  f = (Coor){0.5,-1.0, 0};
  const Coor df = (Coor){1.0, 4.0, 1};
  iCoor P,Q;
  ierr = GridGetBounds(g, &P, &Q); CHKERRQ(ierr);
  CoorToIndex(f, df, g->aabb.lo, &P );
  CoorToIndex(f, df, g->aabb.hi, &Q );
  iCoor S = {0,0,0};
  Coor s;
  for (S.y = P.y; S.y <= Q.y; ++S.y ) {
    for (S.x = P.x; S.x <= Q.x; ++S.x ) {
      IndexToCoor( f, df, S, &s);
      printf("%f, %f; %d, %d\n", s.x, s.y, S.x, S.y );

      findroots(g, s, df, roots);
    }
  }

  ierr = ArrayWrite(roots, 0); CHKERRQ(ierr);
  const int len = ArrayLength( roots );
  printf("len: %d\n", len );

  ierr = ArrayDestroy( roots ); CHKERRQ(ierr);
  ierr = GridDestroy( g ); CHKERRQ(ierr);
  ierr = DCellFinalize(); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode findroots( Grid g, const Coor p0, const Coor df, Array roots )
{
  // upper case: index coor system (P,O,L)
  // lower case: world coor system (p,o)
  const Coor dg = g->d;
  const double TOL = 1e-6;
  const Coor dr = (Coor) { PetscMin( dg.x, df.x), PetscMin( dg.y, df.y), PetscMin( dg.z, df.z) };
  iCoor O;
  Coor P;
  Coor p;
  Coor sol;
  Coor *root;
  PetscErrorCode ierr;

  if( g->is2D ) {
    const PetscReal **phi;
    ierr = GridGet(g, &phi); CHKERRQ(ierr);
    for ( p = p0; p.y < p0.y + df.y; p.y += dr.y ) {
      CoorToIndex2( g->aabb.lo, dg, p, &O, &P);
      if( !GridIndexInBox( g, O) ) continue;
      const Coor L = { P.x - O.x, P.y - O.y, P.z - O.z };
      sol.x = -phi[O.y][O.x] + L.y*phi[O.y][O.x] + phi[O.y][1 + O.x] -
               L.y*phi[O.y][1 + O.x] - L.y*phi[1 + O.y][O.x] + L.y*phi[1 + O.y][1 + O.x];
      if( PetscAbs( sol.x ) < TOL ) continue;
      sol.x = ( -phi[O.y][O.x] + L.y*phi[O.y][O.x] - L.y*phi[1 + O.y][O.x] ) / sol.x;
      if( sol.x < 0 && 1 < sol.x ) continue;
      sol.x *= dg.x;
      if( sol.x < p0.x && p0.x + df.x < sol.x ) continue;
      ierr = ArrayAppend(roots, &root); CHKERRQ(ierr);
      root->x = p.x + sol.x ;
      root->y = p.y;
      root->z = p.z;
    }
    for ( p = p0; p.x < p0.x + df.x; p.x += dr.x ) {
      CoorToIndex2( g->aabb.lo, dg, p, &O, &P);
      if( !GridIndexInBox( g, O) ) continue;
      const Coor L = { P.x - O.x, P.y - O.y, P.z - O.z };
      sol.y = -phi[O.y][O.x] + L.x*phi[O.y][O.x] - L.x*phi[O.y][1 + O.x] +
               phi[1 + O.y][O.x] - L.x*phi[1 + O.y][O.x] + L.x*phi[1 + O.y][1 + O.x];
      if( PetscAbs( sol.y ) < TOL ) continue;
      sol.y = ( -phi[O.y][O.x] + L.x*phi[O.y][O.x] - L.x*phi[O.y][1 + O.x] ) / sol.y;
      if( sol.y < 0 && 1 < sol.y ) continue;
      sol.y *= dg.y;
      if( sol.y < p0.y && p0.y + df.x < sol.y ) continue;
      ierr = ArrayAppend(roots, &root); CHKERRQ(ierr);
      root->x = p.x;
      root->y = p.y + sol.y;
      root->z = p.z;
    }
  } else {

  }
  return 0;
}
