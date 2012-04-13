#include "FluidField.h"

#undef __FUNCT__
#define __FUNCT__ "FluidFieldIntegrateStrainRate"
PetscErrorCode FluidFieldIntegrateStrainRate( DM daV, Vec vecV, DM daE, Vec vecE, PetscReal dh, PetscReal dt )
{
  int i,j,k;
  int xs,ys,zs,
      xm,ym,zm;
  PetscReal dx, dy, dz;
  VelocityVector ***vel;
  StrainTensor ***e;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  dx=dy=dz=dh;

  // Update strain rate
  ierr = DMDAGetCorners(daE,&xs,&ys,&zs,&xm,&ym,&zm); CHKERRQ(ierr);
  ierr = DMDAVecGetArray(daE,vecE,&e); CHKERRQ(ierr);
  ierr = DMDAVecGetArray(daV,vecV,&vel); CHKERRQ(ierr);
  for (k = zs; k < zm; ++k) {
    for (j = ys; j < ym; ++j) {
      for (i = xs; i < xm; ++i) {
        e[k][j][i].xx += dt * ( vel[k][j][i].u - vel[k][j][i].u ) / dx;
        e[k][j][i].yy += dt * ( vel[k][j][i].v - vel[k][j][i].v ) / dy;
        e[k][j][i].zz += dt * ( vel[k][j][i].w - vel[k][j][i].w ) / dz;
        e[k][j][i].xy += dt * 0.5 * ( ( vel[k][j][i].u - vel[k][j][i].u ) / dy +
                                      ( vel[k][j][i].v - vel[k][j][i].v ) / dx );
        e[k][j][i].xz += dt * 0.5 * ( ( vel[k][j][i].u - vel[k][j][i].u ) / dz +
                                      ( vel[k][j][i].w - vel[k][j][i].w ) / dx );
        e[k][j][i].yz += dt * 0.5 * ( ( vel[k][j][i].w - vel[k][j][i].w ) / dy +
                                      ( vel[k][j][i].v - vel[k][j][i].v ) / dz );
      }
    }
  }
  ierr = DMDAVecRestoreArray(daE,vecE,&e); CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(daV,vecV,&vel); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "FluidFieldElasticDivergence"
PetscErrorCode FluidFieldElasticDivergence( DM daE, Vec et, DM daV, Vec rhs, PetscReal dh )
{
  int i,j,k;
  int xs,ys,zs,
      xm,ym,zm;
  PetscReal dx, dy, dz;
  VelocityVector ***b;
  StrainTensor ***e;
  Vec etl; //Local strain tensor
  PetscErrorCode ierr;

  PetscFunctionBegin;
  dx=dy=dz=dh;
  ierr = DMGetLocalVector(daE, &etl); CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(daE,et,INSERT_VALUES,etl); CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(daE,et,INSERT_VALUES,etl); CHKERRQ(ierr);
  ierr = DMDAVecGetArray(daE,et,&e); CHKERRQ(ierr);
  ierr = DMDAVecGetArray(daV,rhs,&b); CHKERRQ(ierr);
  ierr = DMDAGetCorners(daV,&xs,&ys,&zs,&xm,&ym,&zm); CHKERRQ(ierr);
  for (k = zs; k < zm; ++k) {
    for (j = ys; j < ym; ++j) {
      for (i = xs; i < xm; ++i) {
        // Del.e
        b[k][j][i].u += ( e[k][j][i].xx - e[k][j][i].xx ) / dx +
                        ( e[k][j][i].xy - e[k][j][i].xy ) / dy +
                        ( e[k][j][i].xz - e[k][j][i].xz ) / dz;
        b[k][j][i].v += ( e[k][j][i].xy - e[k][j][i].xy ) / dx +
                        ( e[k][j][i].yy - e[k][j][i].yy ) / dy +
                        ( e[k][j][i].yz - e[k][j][i].yz ) / dz;
        b[k][j][i].w += ( e[k][j][i].xz - e[k][j][i].xz ) / dx +
                        ( e[k][j][i].yz - e[k][j][i].yz ) / dy +
                        ( e[k][j][i].zz - e[k][j][i].zz ) / dz;
      }
    }
  }
  ierr = DMDAVecRestoreArray(daV,rhs,&b); CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(daE,et,&e); CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(daE, &etl); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}
