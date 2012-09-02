#include "ImmersedInterfaceMethod.h"
#include "IIM_private.h"

PetscErrorCode IIMIrregularNodes( IIM iim, LevelSet ls );
//#if 0
#undef __FUNCT__
#define __FUNCT__ "IIMIrregularNodes"
PetscErrorCode IIMIrregularNodes( IIM iim, LevelSet ls )
{
  int i, j, k;
  iCoor P,Q;
  Coor pp;
  PetscReal *p = (PetscReal*)&pp;
  PetscReal p0[3], p1[3];
  const Coor   f = iim->f;
  const Coor  df = iim->df;
  const Grid phi = ls->phi;
  const int dims = phi->is2D ? 2 : 3;
  const int faces[3][2] = { {1,2}, {0,2}, {0,1} };
  PetscErrorCode ierr;

  PetscFunctionBegin;
  CoorToIndex( f, df, phi->aabb.lo, &P );
  IndexToCoor( f, df, P, (Coor*)&p0 );
  CoorToIndex( f, df, phi->aabb.hi, &Q );
  Q.x++;Q.y++;Q.z++;
  IndexToCoor( f, df, Q, (Coor*)&p1 );

  for (k = 0; k < dims; ++k) {
    i = faces[k][1];
    j = faces[k][2];
    p[0] = p0[0];
    p[1] = p0[1];
    p[2] = p0[2];
    for ( p[j] = p0[j]; p[j] < p1[j]; p[j] += CoorGet( df, j) ) {
      for ( p[i] = p0[i]; p[i] < p1[i]; p[i] += CoorGet( df, i) ) {
        ierr = GridIntersections( phi, pp, k, iim->roots ); CHKERRQ(ierr);
      }
    }
  }

//  ierr = SpatialIndexInsertPoint(iim->sidx, pt, ir ); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
//#endif

#undef __FUNCT__
#define __FUNCT__ "IIMIrregularNodes2"
PetscErrorCode IIMIrregularNodes2( IIM iim, LevelSet ls, GA ga )
{
  int i,j,k,u,a;
  const int dims = ls->phi->is2D ? 2 : 3;
  const Coor d = iim->df; // Fluid spacing
  Coor eps = ls->phi->d; // LS spacing
  PetscReal **phi;
  Array nodeArray;
  Coor lo, hi;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscLogEventBegin(EVENT_IIMIrregularNodes,0,0,0,0); CHKERRQ(ierr);
  ierr = GridGet(ls->phi, &phi); CHKERRQ(ierr);

  AABB box;
  iCoor p, q;
  eps.x /= 32;
  eps.y /= 32;
  eps.z /= 32;
  const Coor ob[3] = {
      {  1.0, eps.y, eps.z},
      {eps.x,   1.0, eps.z},
      {eps.x, eps.y,   1.0}
  };

  ierr = GridGetBounds(ls->phi, &p, &q); CHKERRQ(ierr);
  for (k = p.z; k < q.z; ++k) {
    for (j = p.y; j < q.y; ++j) {
      for (i = p.x; i < q.x; ++i) {
        Coor lo = {i,j,k};
        for (u = 0; u < dims; ++u) { // face {u, v, w}
          Coor hi = lo;
//          hi[u] += 0.5;

        }
      }
    }
  }
  for (k = p.z; k < q.z; ++k) {
    for (j = p.y; j < q.y; ++j) {
      for (i = p.x; i < q.x; ++i) {
        for (u = 0; u < dims; ++u) { // face {u, v, w}
          for (a = 0; a < dims; ++a) { // direction {x,y,z}
            box.lo.x = d.x * ( i - Tensor1[u][0] - ob[a].x );
            box.hi.x = d.x * ( i - Tensor1[u][0] + ob[a].x );
            box.lo.y = d.y * ( j - Tensor1[u][1] - ob[a].y );
            box.hi.y = d.y * ( j - Tensor1[u][1] + ob[a].y );
            box.lo.z = d.z * ( k - Tensor1[u][2] - ob[a].z );
            box.hi.z = d.z * ( k - Tensor1[u][2] + ob[a].z );

            ierr = SpatialIndexQueryPointsBox( iim->sidx, box, &nodeArray); CHKERRQ(ierr);


          }
        }
      }
    }
  }
  ierr = PetscLogEventEnd(EVENT_IIMIrregularNodes,0,0,0,0); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
