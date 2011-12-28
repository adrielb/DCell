#include "FluidField.h"

#undef __FUNCT__
#define __FUNCT__ "FluidFieldIntegrateStrainRate"
PetscErrorCode FluidFieldIntegrateStrainRate( DA daV, Vec vecV, DA daE, Vec vecE, PetscReal dh, PetscReal dt )
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
  ierr = DAGetCorners(daE,&xs,&ys,&zs,&xm,&ym,&zm); CHKERRQ(ierr);
  ierr = DAVecGetArray(daE,vecE,&e); CHKERRQ(ierr);
  ierr = DAVecGetArray(daV,vecV,&vel); CHKERRQ(ierr);
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
  ierr = DAVecRestoreArray(daE,vecE,&e); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(daV,vecV,&vel); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "FluidFieldElasticDivergence"
PetscErrorCode FluidFieldElasticDivergence( DA daE, Vec et, DA daV, Vec rhs, PetscReal dh )
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
  ierr = DAGetLocalVector(daE, &etl); CHKERRQ(ierr);
  ierr = DAGlobalToLocalBegin(daE,et,INSERT_VALUES,etl); CHKERRQ(ierr);
  ierr = DAGlobalToLocalEnd(daE,et,INSERT_VALUES,etl); CHKERRQ(ierr);
  ierr = DAVecGetArray(daE,et,&e); CHKERRQ(ierr);
  ierr = DAVecGetArray(daV,rhs,&b); CHKERRQ(ierr);
  ierr = DAGetCorners(daV,&xs,&ys,&zs,&xm,&ym,&zm); CHKERRQ(ierr);
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
  ierr = DAVecRestoreArray(daV,rhs,&b); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(daE,et,&e); CHKERRQ(ierr);
  ierr = DARestoreLocalVector(daE, &etl); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}
