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

  iCoor P,Q;
  ierr = GridGetBounds(g, &P, &Q); CHKERRQ(ierr);
  const Coor  f = (Coor){0.5,-1.0,-2.1};
  const Coor df = (Coor){0.5, 3.0,   2};

  CoorToIndex(f, df, g->aabb.lo, &P );
  CoorToIndex(f, df, g->aabb.hi, &Q );
  iCoor S = P;
  Coor s, t;
  for (S.y = P.y; S.y < Q.y; ++S.y ) {
    for (S.x = P.x; S.x < Q.x; ++S.x ) {
      s.x = f.x + S.x * df.x;
      s.y = f.y + S.y * df.y;
      printf("%f, %f; %d, %d\n", s.x, s.y, S.x, S.y );

      t = s;
      for (t.y = s.y; t.y < s.y + df.y; t.y += dg.y ) {
        Coor *root;
        ierr = ArrayAppend(roots, &root); CHKERRQ(ierr);
        root->x = t.x;
        root->y = t.y;
        root->z = t.z;
//        findroots(g, t, df, roots);
      }

      t = s;
      for (t.x = s.x; t.x < s.x + df.x; t.x += dg.x ) {
//        findroots(g, t, df, roots);
      }
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

PetscErrorCode findroots( Grid g, const Coor p, const Coor dr, Array roots )
{
  if( !AABBPointInBox(g->aabb, p) )
    return 0;
  // upper case: index coor system (P,O,L)
  // lower case: world coor system (p,o)
  const Coor s = g->aabb.lo;
  const Coor d = g->d;
  const Coor P = {
      (p.x - s.x) / d.x,
      (p.y - s.y) / d.y,
      (p.z - s.z) / d.z };
  const iCoor O = { (int)floor( P.x ), (int)floor( P.y ), (int)floor( P.z ) };
  const  Coor LL = { P.x - O.x, P.y - O.y, P.z - O.z };
  const double TOL = 1e-6;
  Coor L;
  Coor sol;
  Coor *root;
  PetscErrorCode ierr;

  if( g->is2D ) {
    const PetscReal **phi;
    ierr = GridGet(g, &phi); CHKERRQ(ierr);

    L = LL;
    for (; L.y < 1.0; L.y += dr.y ) {
      sol.x = -phi[O.y][O.x] + L.y*phi[O.y][O.x] + phi[O.y][1 + O.x] -
               L.y*phi[O.y][1 + O.x] - L.y*phi[1 + O.y][O.x] + L.y*phi[1 + O.y][1 + O.x];
      if( PetscAbs( sol.x ) > TOL ) {
        sol.x = ( -phi[O.y][O.x] + L.y*phi[O.y][O.x] - L.y*phi[1 + O.y][O.x] ) / sol.x;
        if( 0 <= sol.x && sol.x < 1 ) {
          ierr = ArrayAppend(roots, &root); CHKERRQ(ierr);
          root->x = O.x + sol.x ;
          root->y = O.y + L.y;
          root->z = O.z + L.z;
        }
      }
    }

    L = LL;
    for (; L.x < 1.0; L.x += dr.x ) {
      sol.y = -phi[O.y][O.x] + L.x*phi[O.y][O.x] - L.x*phi[O.y][1 + O.x] +
               phi[1 + O.y][O.x] - L.x*phi[1 + O.y][O.x] + L.x*phi[1 + O.y][1 + O.x];
      if( PetscAbs( sol.y ) > TOL ) {
        sol.y = ( -phi[O.y][O.x] + L.x*phi[O.y][O.x] - L.x*phi[O.y][1 + O.x] ) / sol.y;
        if( 0 <= sol.y && sol.y < 1 ) {
          ierr = ArrayAppend(roots, &root); CHKERRQ(ierr);
          root->x = O.x + L.x;
          root->y = O.y + sol.y;
          root->z = O.z + L.z;
        }
      }
    }
  } else {

  }
  return 0;
}
