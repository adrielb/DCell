#include "LevelSetMethod.h"

PetscErrorCode SetVelocity(Grid velgrid);

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc, char **args)
{
  PetscErrorCode  ierr;
  PetscFunctionBegin;
  ierr = PetscInitialize(&argc, &args, (char *) 0, ""); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "Start %s\n", __FILE__); CHKERRQ(ierr);

  int d1 = 16;
  PetscReal dx = 2./(d1-1);
  Coor dh = {dx, dx, dx};

  LevelSet ls;
  Coor pos = {0,0,0};
//  ierr = LevelSetInitializeToStar2D(dh, pos, 0.2, 0.5, 5, &ls); CHKERRQ(ierr);
  ierr = LevelSetInitializeToCircle(dh, pos, 0.8, &ls); CHKERRQ(ierr);
  ierr = LevelSetInitializeParticles(ls); CHKERRQ(ierr);

  Grid velgrid;
  Grid phi = ls->phi;
  iCoor size = phi->n;
  printf("ls grid: [%d,%d,%d]\n", size.x, size.y, size.z);
  size.x *= 10;
  size.y *= 10;
  size.z *= 10;
  printf("velgrid: [%d,%d,%d]\n", size.x, size.y, size.z);
  ierr = GridCreate(phi->d,(iCoor){-size.x/2,-size.y/2,-size.z/2},size,phi->is2D?2:3,&velgrid); CHKERRQ(ierr);
  ierr = GridSetName(velgrid,"velgrid"); CHKERRQ(ierr);

  ierr = SetVelocity(velgrid); CHKERRQ(ierr);
  ierr = GridWrite(velgrid,0); CHKERRQ(ierr);
  ierr = GridWrite(ls->phi,0); CHKERRQ(ierr);
  ierr = LevelSetWriteIrregularNodeList( ls, 0 ); CHKERRQ(ierr);
  ierr = ParticleLSWriteParticles(ls->pls, 0); CHKERRQ(ierr);
  ls->AdvectThres = 1;

  PetscReal dt = dx/5;
  PetscInt t, tmax = 2*PETSC_PI*d1/2.;
  tmax = 18*5;
  for (t = 1; t < tmax; ++t) {
    printf("t: %d\t:\t%d\t%f\n",t,tmax,(100.*t) / tmax);
    ierr = ls->Advect( ls, velgrid, dt); CHKERRQ(ierr);
    ierr = LevelSetWriteIrregularNodeList( ls, t ); CHKERRQ(ierr);
    ierr = GridWrite(ls->phi,t); CHKERRQ(ierr);
    ierr = ParticleLSWriteParticles(ls->pls, t); CHKERRQ(ierr);
  }

  ierr = LevelSetDestroy(ls); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "End %s\n", __FILE__); CHKERRQ(ierr);
  ierr = PetscFinalize(); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode SetVelocity(Grid velgrid)
{
  int i,j;
  PetscReal x,y,r;
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
      r = PetscSqrtScalar(x*x+y*y);
      r = r == 0.0 ? 1 : r;
      vel[j][i][0] = x / r;
//      vel[j][i][0] = -0.1*PetscSign(x);

      x = dx * i;
      y = dx * (j - Tensor1[0][1]);
      r = PetscSqrtScalar(x*x+y*y);
      r = r == 0.0 ? 1 : r;
      vel[j][i][1] = y / r;
//      vel[j][i][1] = -0.1*PetscSign(y);
    }
  }
  return 0;
}

/*
  DA da;
  ierr = DACreate2d(PETSC_COMM_WORLD,DA_NONPERIODIC,DA_STENCIL_BOX,d1,d1,1,1,2,1,0,0,&da); CHKERRQ(ierr);
  Vec vecVel;
  ierr = DACreateGlobalVector(da,&vecVel); CHKERRQ(ierr);
  PetscReal ***vel;
  ierr = DAVecGetArrayDOF(da,vecVel,&vel); CHKERRQ(ierr);
  ierr = DAVecRestoreArrayDOF(da,vecVel,&vel); CHKERRQ(ierr);
      X = x * (E.x - S.x) / (d1 - 1) + S.x;
      Y = y * (E.y - S.y) / (d1 - 1) + S.y;
 */
