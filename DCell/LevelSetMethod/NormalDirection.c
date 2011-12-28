#include "LevelSetMethod.h"

#undef __FUNCT__
#define __FUNCT__ "LevelSetNormalDirection"
PetscErrorCode LevelSetNormalDirection( LevelSet ls, Coor X, Coor *n )
{
  const Coor dh = ls->phi->d;
  PetscReal h;

  PetscErrorCode ierr;
  // TODO: profile. if test vs function pointer for 2D/3D
  if( ls->phi->is2D ) {
    PetscReal **phi;
    ierr = GridGet(ls->phi,&phi); CHKERRQ(ierr);
    n->x = Bilinear2D(GridFunction2D_DerivX, phi, dh, X.x, X.y);
    n->y = Bilinear2D(GridFunction2D_DerivY, phi, dh, X.x, X.y);
    h = sqrt( n->x*n->x + n->y*n->y );
    n->x = n->x / h;
    n->y = n->y / h;
    n->z = 0;
  } else {
    PetscReal ***phi;
    ierr = GridGet(ls->phi,&phi); CHKERRQ(ierr);
    n->x = Bilinear3D(GridFunction3D_DerivX, phi, dh, X);
    n->y = Bilinear3D(GridFunction3D_DerivY, phi, dh, X);
    n->z = Bilinear3D(GridFunction3D_DerivZ, phi, dh, X);
    h = sqrt( n->x*n->x + n->y*n->y + n->z*n->z );
    n->x = n->x / h;
    n->y = n->y / h;
    n->z = n->z / h;
  }

  PetscFunctionReturn(0);
}

inline PetscReal LevelSetDiracDelta2D( PetscReal **phi, const Coor dh, const Coor X )
{
  const PetscReal eps = 1.5;
  const PetscReal phival = Bilinear2D(GridFunction2D_Identity,phi,dh,X.x,X.y);
  if( phival < -eps || phival > eps ) return 0;
  const PetscReal px = Bilinear2D(GridFunction2D_DerivX,phi,dh,X.x,X.y);
  const PetscReal py = Bilinear2D(GridFunction2D_DerivY,phi,dh,X.x,X.y);
  const PetscReal gradmag = PetscSqrtScalar( px*px + py*py );
  return gradmag * ( 1 + cos( PETSC_PI * phival / eps) ) / (2*eps);
}
