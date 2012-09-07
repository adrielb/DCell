#include "LevelSetMethod.h"

PetscErrorCode LevelSetNormalDirection( LevelSet ls, Coor X, Coor *n )
{
  Grid phi = ls->phi;
  PetscReal h;

  // TODO: profile. if test vs function pointer for 2D/3D
  if( ls->phi->is2D ) {
    n->x = GridBilinear( phi, GridFunction2D_DerivX, X );
    n->y = GridBilinear( phi, GridFunction2D_DerivY, X );
    h = sqrt( n->x*n->x + n->y*n->y );
    n->x = n->x / h;
    n->y = n->y / h;
    n->z = 0;
  } else {
    n->x = GridBilinear( phi, GridFunction3D_DerivX, X );
    n->y = GridBilinear( phi, GridFunction3D_DerivY, X );
    n->z = GridBilinear( phi, GridFunction3D_DerivZ, X );
    h = sqrt( n->x*n->x + n->y*n->y + n->z*n->z );
    n->x = n->x / h;
    n->y = n->y / h;
    n->z = n->z / h;
  }

  return 0;
}

inline PetscReal LevelSetDiracDelta2D( Grid phi, const Coor X )
{
  const PetscReal eps = 1.5;
  PetscReal phival;
  GridInterpolate(phi, X, &phival );
  if( phival < -eps || phival > eps ) return 0;
  const PetscReal px = GridBilinear( phi, GridFunction2D_DerivX, X );
  const PetscReal py = GridBilinear( phi, GridFunction2D_DerivY, X );
  const PetscReal gradmag = PetscSqrtScalar( px*px + py*py );
  return gradmag * ( 1 + cos( PETSC_PI * phival / eps) ) / (2*eps);
}

inline PetscReal LevelSetDiracDelta3D( Grid phi, const Coor X )
{
  const PetscReal eps = 1.5;
  PetscReal phival;
  GridInterpolate(phi, X, &phival );
  if( phival < -eps || phival > eps ) return 0;
  const PetscReal px = GridBilinear( phi, GridFunction3D_DerivX, X);
  const PetscReal py = GridBilinear( phi, GridFunction3D_DerivY, X);
  const PetscReal pz = GridBilinear( phi, GridFunction3D_DerivZ, X);
  const PetscReal gradmag = PetscSqrtScalar( px*px + py*py + pz*pz );
  return gradmag * ( 1 + cos( PETSC_PI * phival / eps) ) / (2*eps);
}
