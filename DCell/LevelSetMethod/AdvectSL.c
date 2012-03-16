#include "LevelSetMethod.h"
#include "LSM_private.h"

#undef  __FUNCT__
#define __FUNCT__ "LevelSetAdvectSL"
PetscErrorCode LevelSetAdvectSL(LevelSet ls, Grid velgrid, PetscReal dt)
{
  int b;
  iCoor *band;
  Coor X; // Grid point
  Coor S; // Projected grid point
  Coor V; // Interpolated velocity
  Coor dh = ls->phi->d;
  PetscReal ***vel, **phi, **phi0;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscLogEventBegin(EVENT_LevelSetAdvectSL,0,0,0,0); CHKERRQ(ierr);
  ierr = GridGet(velgrid,&vel); CHKERRQ(ierr);
  ierr = GridGet(ls->phi,&phi); CHKERRQ(ierr);
  ierr = GridGet(ls->phi0,&phi0); CHKERRQ(ierr);

  for ( b = 0; b < ArrayLength(ls->band); ++b) {
    ierr = ArrayGet(ls->band,b,&band); CHKERRQ(ierr);
    X.x = band->x;
    X.y = band->y;
    // Interpolate velocity at grid point
    ierr = InterpolateVelocity2D( 0, vel, X, &V ); CHKERRQ(ierr);
    // Project grid point back in time
    S.x = X.x - V.x * dt / dh.x;
    S.y = X.y - V.y * dt / dh.y;
    // Interpolate phi at projection
    phi[band->y][band->x] = Bilinear2D(GridFunction2D_Identity,phi0,dh, S.x, S.y );
  } // for b in band
  ierr = PetscLogEventEnd(EVENT_LevelSetAdvectSL,0,0,0,0); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef  __FUNCT__
#define __FUNCT__ "LevelSetAdvectSL_3D"
PetscErrorCode LevelSetAdvectSL_3D(LevelSet ls, Grid velgrid, PetscReal dt)
{
  int b;
  int len = ArrayLength(ls->band);
  iCoor *band = ArrayGetData(ls->band);
  Coor X; // Grid point
  Coor S; // Projected grid point
  Coor V; // Interpolated velocity
  Coor dh = ls->phi->d;
  PetscReal ****vel, ***phi, ***phi0;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscLogEventBegin(EVENT_LevelSetAdvectSL,0,0,0,0); CHKERRQ(ierr);
  ierr = GridGet(velgrid,&vel); CHKERRQ(ierr);
  ierr = GridGet(ls->phi ,&phi); CHKERRQ(ierr);
  ierr = GridGet(ls->phi0,&phi0); CHKERRQ(ierr);

  for ( b = 0; b < len; ++b) {
    X.x = band[b].x;
    X.y = band[b].y;
    X.z = band[b].z;

    ierr = InterpolateVelocity3D( 0, vel, X, &V ); CHKERRQ(ierr);

    S.x = X.x - V.x * dt / dh.x;
    S.y = X.y - V.y * dt / dh.y;
    S.z = X.z - V.z * dt / dh.z;

    phi[band[b].z][band[b].y][band[b].x] = Bilinear3D(GridFunction3D_Identity, phi0, dh, S );
  } // for b in band
  ierr = PetscLogEventEnd(EVENT_LevelSetAdvectSL,0,0,0,0); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
