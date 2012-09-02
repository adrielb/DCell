#include "ImmersedInterfaceMethod.h"
#include "IIM_private.h"

PetscErrorCode findroots( Grid g, const Coor p, const Coor dr, Array roots );
PetscErrorCode myroot( Grid g, Coor p, int axis, Array roots );

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc, char **args)
{
  PetscErrorCode  ierr;
  PetscFunctionBegin;
  ierr = DCellInit(); CHKERRQ(ierr);

  const Coor dg = {0.5,0.5,0.5};
  const iCoor size = {3,3,0};
  const iCoor pos = {0,0,0};
  Grid g;
  ierr = GridCreate( dg, pos, size, 1, &g); CHKERRQ(ierr);
  PetscReal **val;
  ierr = GridGet(g, &val); CHKERRQ(ierr);
  val[2][0] = 1; val[2][1] = 1; val[2][2] = 15;
  val[1][0] = 1; val[1][1] =  -1; val[1][2] = 1;
  val[0][0] = 0.05; val[0][1] = 0.5; val[0][2] = 1;
  ierr = GridWrite(g, 0); CHKERRQ(ierr);

  IIM iim;
  ierr = IIMCreate( g->is2D, &iim); CHKERRQ(ierr);

  const Coor  f = (Coor){-0.1,  -0.1, 0};
  const Coor df = (Coor){ 0.06, 0.06, 1};

  ierr = IIMSetFluidCoordinateSystem( iim, f, df); CHKERRQ(ierr);

  iCoor P,Q;
  CoorToIndex(f, df, g->aabb.lo, &P );
  Coor pp;
  PetscReal *p= (PetscReal*)&pp, p0[3], p1[3];
  IndexToCoor( f, df, P, (Coor*)&p0 );
  CoorToIndex(f, df, g->aabb.hi, &Q );
  Q.x++;Q.y++;Q.z++;
  IndexToCoor( f, df, Q, (Coor*)&p1 );
  printf( "p0: %f, %f, %f\n", p0[0], p0[1], p0[2] );
  printf( "p1: %f, %f, %f\n", p1[0], p1[1], p1[2] );
  int i, j, k;
//  const int dims = g->is2D ? 2 : 3;
  const int dims = 3;
  const int faces[3][3] = { {0,1,2}, {1,0,2}, {2,0,1} };
  for (k = 0; k < dims; ++k) {
    i = faces[k][1];
    j = faces[k][2];
    printf("========================================\n" );
    printf("%d %d\n", i, j);
    p[0] = p0[0];
    p[1] = p0[1];
    p[2] = p0[2];
    for ( p[j] = p0[j]; p[j] < p1[j]; p[j] += CoorGet( df, j) ) {
      for ( p[i] = p0[i]; p[i] < p1[i]; p[i] += CoorGet( df, i) ) {
//        printf( "{%f, %f, %f},\n", p[0], p[1], p[2] );
        ierr = GridIntersections( g, pp, k, iim->roots ); CHKERRQ(ierr);
      }
    }
  }

  ierr = ArrayWrite( iim->roots, 0); CHKERRQ(ierr);
  const int len = ArrayLength( iim->roots );
  printf("len: %d\n", len );

  ierr = IIMDestroy( iim ); CHKERRQ(ierr);
  ierr = GridDestroy( g ); CHKERRQ(ierr);
  ierr = DCellFinalize(); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


