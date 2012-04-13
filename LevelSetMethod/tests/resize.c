#include "LevelSetMethod.h"

PetscErrorCode CreateVelocityGrid( LevelSet ls, Grid velgrid, int sign);
PetscErrorCode LevelSetAdvectSL( LevelSet ls, Grid velgrid, PetscReal dt);

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc, char **args)
{
  PetscErrorCode  ierr;
  PetscFunctionBegin;
  ierr = DCellInit(); CHKERRQ(ierr);

  int d1 = 32;
  Coor dh = {1./(d1-1), 1./(d1-1), 1./(d1-1)};

  PetscReal radius = 0.4;
  LevelSet ls;
  ierr = LevelSetInitializeToStar2D( dh, (Coor){0,0,0}, radius, 0.4, 8, &ls); CHKERRQ(ierr);
  ierr = LevelSetWriteIrregularNodeList(ls->irregularNodes,0); CHKERRQ(ierr);

  Grid velgrid;
  ierr = GridCreate(dh,ls->phi->p,ls->phi->n,2,&velgrid); CHKERRQ(ierr);

  PetscReal maxVel = sqrt(2);
  PetscReal tend = 4*radius/maxVel;
  PetscReal CFL = .1;
  PetscReal dt = CFL * dh.x / maxVel;
  ierr = GridWrite(ls->phi,0); CHKERRQ(ierr);
  ierr = ArrayWrite(ls->band,"band",0); CHKERRQ(ierr);
  int sign = 1;
  for (int t = 1; t*dt < tend; ++t) {
    if( t*dt > tend/2 ) sign = -1;
    ierr = LevelSetWriteIrregularNodeList(ls->irregularNodes,t); CHKERRQ(ierr);
    ierr = CreateVelocityGrid(ls, velgrid, sign); CHKERRQ(ierr);
    ierr = LevelSetCFLIncrement( ls, velgrid, dt ); CHKERRQ(ierr);
    ierr = LevelSetAdvectSL( ls, velgrid, dt ); CHKERRQ(ierr);
    ierr = LevelSetUpdateIrregularNodeList( ls ); CHKERRQ(ierr);
    if( ls->CFLcount > ls->CFLthres ) {
      ls->CFLcount = 0;
      ierr = LevelSetReinitialize(ls); CHKERRQ(ierr);
    }
    ierr = GridWrite(ls->phi,t); CHKERRQ(ierr);
    ierr = ArrayWrite(ls->band,"band",t); CHKERRQ(ierr);
    printf("t: %1.3f %1.3f \t %d %1.1f\t\t\t%f\n",t*dt, tend, t, tend/dt, ls->CFLcount);
  }

  ierr = LevelSetDestroy(ls); CHKERRQ(ierr);
  ierr = DCellFinalize(); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode CreateVelocityGrid( LevelSet ls, Grid velgrid, int sign)
{
  PetscErrorCode ierr;
  ierr = GridResize(velgrid,ls->phi->p,ls->phi->n); CHKERRQ(ierr);
  int i,j;
  iCoor p,q;
  PetscReal ***vel;
  ierr = GridGetBounds(velgrid,&p,&q); CHKERRQ(ierr);
  ierr = GridGet(velgrid,&vel); CHKERRQ(ierr);
  for (j = p.y; j < q.y; ++j) {
    for (i = p.x; i < q.x; ++i) {
      vel[j][i][0] = sign;
      vel[j][i][1] = sign;
    }
  }
  return 0;
}
