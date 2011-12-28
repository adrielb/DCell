#include "LevelSetMethod.h"

PetscErrorCode SetVelocity(Grid velgrid, PetscReal t);

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc, char **args) {

  PetscErrorCode ierr;
  ierr = DCellInit(); CHKERRQ(ierr);
  PetscReal dx = 1;
  PetscReal radius = 10;
  Coor center = (Coor){0,0,0};
  Coor dh = {dx,dx,0};
  LevelSet ls;
  ierr = LevelSetInitializeToStar2D(dh,center,radius,2,5,&ls); CHKERRQ(ierr);
  ierr = LevelSetWriteIrregularNodeList( ls->irregularNodes, 0 ); CHKERRQ(ierr);

  Grid velgrid;
  Grid phi = ls->phi;
  ierr = GridCreate(phi->d,phi->p,phi->n,phi->is2D?2:3,&velgrid); CHKERRQ(ierr);
  ierr = GridSetName(velgrid,"velgrid"); CHKERRQ(ierr);

  PetscReal t = 0;
  PetscReal dt = dx / radius;
  int ti;
  for (t = 0; t < 1; ) {
    printf("ti: %d\n",ti, t);
    ierr = SetVelocity(velgrid, t); CHKERRQ(ierr);
    ierr = LevelSetAdvectSLRK2HalfStep( ls, velgrid, dt/2 ); CHKERRQ(ierr);
    ierr = SetVelocity(velgrid, t+dt/2); CHKERRQ(ierr);
    ierr = LevelSetAdvectSLRK2FullStep( ls, velgrid, dt); CHKERRQ(ierr);
    ti++;
    ierr = LevelSetWriteIrregularNodeList( ls->irregularNodes, ti ); CHKERRQ(ierr);
    t = t + dt;
  }

  ierr = LevelSetDestroy(ls); CHKERRQ(ierr);
	ierr = DCellFinalize(); CHKERRQ(ierr);
  return 0;
}

PetscErrorCode SetVelocity(Grid velgrid, PetscReal t)
{
  int i,j;
  PetscReal x,y;
  iCoor p,q;
  PetscReal ***vel;
  PetscReal dx = velgrid->d.x;
  PetscErrorCode ierr;

  ierr = GridGet(velgrid,&vel); CHKERRQ(ierr);
  ierr = GridGetBounds(velgrid,&p,&q); CHKERRQ(ierr);
  for ( j = p.y; j < q.y; ++j) {
    for ( i = p.x; i < q.x; ++i) {
      x = dx * (i-Tensor1[0][0]);
      y = dx * j;
      vel[j][i][0] = -y;

      x = dx * i;
      y = dx * (j - Tensor1[0][1]);
      vel[j][i][1] = x;
    }
  }
  return 0;
}
