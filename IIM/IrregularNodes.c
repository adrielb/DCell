#include "ImmersedInterfaceMethod.h"
#include "IIM_private.h"

PetscErrorCode IIMIrregularNodes( IIM iim, LevelSet ls, int GA );

#undef __FUNCT__
#define __FUNCT__ "IIMIrregularNodes"
PetscErrorCode IIMIrregularNodes( IIM iim, LevelSet ls, GA ga )
{
  int i,j,k,u,a;
  const int dims = ls->phi->is2D ? 2 : 3;
  const Coor d = iim->dh; // Fluid spacing
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


#undef __FUNCT__
#define __FUNCT__ "GridRoots"
PetscErrorCode GridRoots( Grid g, Coor s, Coor dF, Array roots )
{
  int i,j,k;
  const Coor dg = g->d;
  const Coor dx = (Coor){ dg.x < dF.x ? dg.x : dF.x,
                          dg.y < dF.y ? dg.y : dF.y,
                          dg.z < dF.z ? dg.z : dF.z };

  PetscErrorCode ierr;

  PetscFunctionBegin;
  for (j = 0; j < 10; ++j) {
    for (i = 0; i < 10; ++i) {

    }
  }
  PetscFunctionReturn(0);
}

