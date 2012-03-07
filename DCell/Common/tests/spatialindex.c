#include "Common.h"

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc, char **args)
{
  PetscErrorCode  ierr;
  PetscFunctionBegin;
  ierr = DCellInit(); CHKERRQ(ierr);
  Coor lo = {0,0,0};
  Coor hi = {4,4,4};
  PetscReal dx = 0.02;
  Coor dh = {dx,dx,dx};
  SpatialIndex sidx;
  ierr = SpatialIndexCreate(lo, hi, dh,  &sidx); CHKERRQ(ierr);

  PetscRandom rnd;
  ierr = PetscRandomCreate(PETSC_COMM_WORLD,&rnd); CHKERRQ(ierr);
  ierr = PetscRandomSetInterval(rnd,lo.x-2,hi.x+2); CHKERRQ(ierr);
//  ierr = PetscRandomSetInterval(rnd,1.2,1.5); CHKERRQ(ierr);
  ierr = PetscRandomSetType(rnd,PETSCRAND48); CHKERRQ(ierr);

  int n = 64*1000*1000;
  Coor *pt;
  int i, j, k;
  Coor center = {1.1, 1.1, 1.1};
  PetscReal r = 0.08;
  const int Np = 10;
  int len;
  Coor *list[Np];

  ierr = PetscMalloc(sizeof(Coor)*n,&pt); CHKERRQ(ierr);
  for (i = 0; i < n; ++i) {
    ierr = PetscRandomGetValueReal(rnd,&pt[i].x); CHKERRQ(ierr);
    ierr = PetscRandomGetValueReal(rnd,&pt[i].y); CHKERRQ(ierr);
    ierr = PetscRandomGetValueReal(rnd,&pt[i].z); CHKERRQ(ierr);
    ierr = SpatialIndexInsertPoint(sidx,pt[i],&pt[i]); CHKERRQ(ierr);
  }

  // loop to test memory leaks
  for (j = 0; j < 100; ++j) {
    printf("j: %d\n",j);
    for (k = 0; k < 100; ++k) {
      ierr = PetscRandomGetValueReal(rnd,&center.x); CHKERRQ(ierr);
      ierr = PetscRandomGetValueReal(rnd,&center.y); CHKERRQ(ierr);
      ierr = PetscRandomGetValueReal(rnd,&center.z); CHKERRQ(ierr);
      ierr = SpatialIndexQueryPoints(sidx, center, r, Np, &len, &list); CHKERRQ(ierr);
    }

    for (i = 0; i < len; ++i) {
      //printf("{%f,%f,%f},", list[i]->x, list[i]->y, list[i]->z);
    }
  }

  ierr = SpatialIndexClear( sidx ); CHKERRQ(ierr);

  ierr = SpatialIndexDestroy( sidx ); CHKERRQ(ierr);
  ierr = DCellFinalize(); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
