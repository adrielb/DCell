#include "ImmersedInterfaceMethod.h"
#include "IIM_private.h"

#undef __FUNCT__
#define __FUNCT__ "IIMGridIntersections"
PetscErrorCode IIMGridIntersections( IIM iim, const Grid g, const Coor p, const int axis )
{
  Coor X0, X1;
  PetscReal *x0 = &((PetscReal*)&X0)[axis];
  PetscReal *x1 = &((PetscReal*)&X1)[axis];
  PetscReal v0, v1, x;
  PetscReal dh = CoorGet( g->d, axis) / 10;
  PetscReal q = CoorGet(g->aabb.hi,axis);
  Coor root;
  IIMIrregularNode *newnode;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  q += dh;
  for ( X0 = p; *x0 < q; *x0 += dh ) {
    X1 = X0;
    *x1 += dh;
    GridInterpolate( g, X0, &v0 );
    GridInterpolate( g, X1, &v1 );
    x = (v0*(*x1) - v1*(*x0)) / (v0 - v1);
    if( x != x || x < *x0 || x >= *x1 ) continue;
    root = X0;
    CoorSet( &root, axis, x);
    if( !AABBPointInBox( g->aabb, root) ) continue;
    ierr = ArrayAppend(iim->irregularNodes, &newnode); CHKERRQ(ierr);
    newnode->X = root;
    ierr = SpatialIndexInsertPoint( iim->sidx, root, newnode ); CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "IIMUpdateIrregularNodes"
PetscErrorCode IIMUpdateIrregularNodes( IIM iim, LevelSet ls )
{
  int i, j;
  int axis;
  iCoor P,Q;
  Coor p;
  Coor p0, p1;
  const Coor   f = iim->f;
  const Coor  df = iim->df;
  const Grid phi = ls->phi;
  const int dims = phi->is2D ? 2 : 3;
  const int faces[3][2] = { {1,2}, {0,2}, {0,1} };
  const AABB aabb = phi->aabb;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = SpatialIndexClear( iim->sidx ); CHKERRQ(ierr);
  ierr = SpatialIndexSetDomain( iim->sidx, aabb.lo, aabb.hi, df ); CHKERRQ(ierr);

  CoorToIndex( f, df, phi->aabb.lo, &P );
  IndexToCoor( f, df, P, &p0 );
  CoorToIndex( f, df, phi->aabb.hi, &Q );
  Q.x++;Q.y++;Q.z++;
  IndexToCoor( f, df, Q, &p1 );

  for (axis = 0; axis < dims; ++axis) {
    i = faces[axis][0];
    j = faces[axis][1];
    p = p0;
    for (   CoorSet(&p,j, CoorGet(p0,j)); CoorGet(p,j) < CoorGet(p1,j); CoorSet(&p,j, CoorGet(p,j)+CoorGet(df,j)) ) {
      for ( CoorSet(&p,i, CoorGet(p0,i)); CoorGet(p,i) < CoorGet(p1,i); CoorSet(&p,i, CoorGet(p,i)+CoorGet(df,i)) ) {
        ierr = IIMGridIntersections( iim, phi, p, axis ); CHKERRQ(ierr);
      }
    }
  }
  PetscFunctionReturn(0);
}

#if 0
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
#endif
