#include "ImmersedInterfaceMethod.h"

inline PetscReal IIMVelocityCorrection( Coor X, IrregularNode *n, Coor *vel )
{
  PetscReal s,a;
  Coor o0;
  Coor o1 = {o0.x+1, o0.y+1, o0.z+1};
  xs = x - n->X.x;
  ys = y - n->X.y;
  zs = z - n->X.z;

  xi = 1 - PetscAbs( xs );
  yi = 1 - PetscAbs( ys );
  zi = 1 - PetscAbs( zs );

  if( xi < 0 ) return 0;
  if( yi < 0 ) return 0;
  if( zi < 0 ) return 0;

  if( n->shift == CELL_CENTER ) {
    s = x[n->axis];
    a = ((PetscReal)n->X)[n->axis];
  }
  const PetscReal C = s < a ? (1-a) * s : a * (1-s);

  vel[dof] += xi * yi * zi * C * uj;
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
