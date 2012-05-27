#include "ImmersedInterfaceMethod.h"

inline void IIMVelocityCorrection( Coor X, IrregularNode *n, Coor *vel )
{
  // orthogonal projection node
  if( n->axis == -1 ) return;

  // stencil intersection
  int uvw = n->shift == CELL_CENTER ? n->axis : n->shift - U_FACE;

  PetscReal s,a;
  Coor O0;
  Coor O1 = {O0.x+1, O0.y+1, O0.z+1};
  PetscReal o0=0, o1=1;
  PetscReal *x = &X.x;
  PetscReal *v = &vel->x;
  PetscReal *nX= &n->X.x;
  PetscInt  *p = &n->pos.x;
  PetscReal xs, ys, zs;
  Coor Xi;
  PetscReal *xi = &Xi.x;


  Tensor1[uvw][0];

  xs = X.x - n->X.x;
  ys = X.y - n->X.y;
  zs = X.z - n->X.z;

  Xi.x = 1 - PetscAbs( xs );
  Xi.y = 1 - PetscAbs( ys );
  Xi.z = 1 - PetscAbs( zs );

  if( Xi.x < 0 ) return;
  if( Xi.y < 0 ) return;
  if( Xi.z < 0 ) return;

  if( n->shift == CELL_CENTER ) {
    s = x[n->axis];
    a = nX[n->axis];
    xi[n->axis] = 1;
    o0 = p[n->axis] - 0.5;
    o1 = p[n->axis] + 0.5;
  } else {

  }

  if( s < o0 ) return;
  if( s > o1 ) return;
  PetscReal C = s < a ? (o1-a) * (s-o0) : (a-o0) * (o1-s);

  PetscReal uj = 1;

  v[uvw] += Xi.x * Xi.y * Xi.z * C * uj;
//  v[uvw] = a;
}

PetscErrorCode IIMCorrectVelocity( IIM iim, const Coor X, Coor *vel )
{
  int i;
  int len;
  const int MAXLEN = 64;
  const PetscReal radius = 2.0;
  IrregularNode *nodes[MAXLEN];
  PetscErrorCode ierr;

  ierr = SpatialIndexQueryPoints( iim->sidx, X, radius, MAXLEN, &len, (void*)&nodes); CHKERRQ(ierr);
  for (i = 0; i < len; ++i) {
    IIMVelocityCorrection( X, nodes[i], vel );
  }
  return 0;
}
