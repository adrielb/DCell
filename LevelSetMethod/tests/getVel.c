#include "LevelSetMethod.h"

PetscErrorCode LevelSetGetVelocity(LevelSet ls, int ga, Grid velgrid);

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc, char **args) {

  PetscErrorCode ierr;
  ierr = DCellInit(); CHKERRQ(ierr);
  iCoor dims = {10,10,0};
  int dof = 3; // [u v p]
  DA daV;
  ierr = DACreate2d(PETSC_COMM_WORLD,DA_NONPERIODIC,DA_STENCIL_STAR,
            dims.x,dims.y, PETSC_DECIDE,PETSC_DECIDE, dof,1, 0,0, &daV); CHKERRQ(ierr);
  int ga;
  ierr = GACreate( daV, &ga); CHKERRQ(ierr);

  Vec global;
  ierr = DACreateGlobalVector(daV, &global); CHKERRQ(ierr);
  PetscReal ***g;
  ierr = DAVecGetArrayDOF(daV,global,&g); CHKERRQ(ierr);
  DALocalInfo info;
  ierr = DAGetLocalInfo(daV,&info); CHKERRQ(ierr);
  for (int j = info.ys; j < info.ys+info.ym; ++j) {
    for (int i = info.xs; i < info.xs+info.xm; ++i) {
      g[j][i][0] = i;
      g[j][i][1] = j;
    }
  }
  ierr = DAVecRestoreArrayDOF(daV,global,&g); CHKERRQ(ierr);
  ierr = GAPutVec(global,ga); CHKERRQ(ierr);

  int rank;
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank ); CHKERRQ(ierr);
  if( rank == 0 )
  {
    LevelSet ls;
    Coor dh = {1,1,0};
    iCoor pos = {-5,-5,0};
    iCoor size = {19,19,0};
    ierr = LevelSetCreate(dh,pos,size,&ls); CHKERRQ(ierr);

    int dof = 2;
    Grid velgrid;
    ierr = GridCreate(dh,pos,size,dof,&velgrid); CHKERRQ(ierr);

    ierr = LevelSetGetVelocity(ls,ga,velgrid); CHKERRQ(ierr);

    for (int i = 0; i < velgrid->SIZE; ++i) {
      printf( "%3.0f ", velgrid->v1[i] );
    }
    printf("\n\n");

    PetscReal ***grid;
    ierr = GridGet(velgrid, &grid); CHKERRQ(ierr);
    iCoor p, q;
    ierr = GridGetBounds(velgrid, &p, &q); CHKERRQ(ierr);
    for (int j = p.y; j < q.y; ++j) {
      for (int i = p.x; i < q.x; ++i) {
        printf("%3.f\t%3.f\t:", grid[j][i][0], grid[j][i][1]);
      }
      printf("\n");
    }

    ierr = GridDestroy(velgrid); CHKERRQ(ierr);
    ierr = LevelSetDestroy(ls); CHKERRQ(ierr);
  }

  GA_Destroy(ga);
  ierr = DADestroy(daV); CHKERRQ(ierr);
	ierr = DCellFinalize(); CHKERRQ(ierr);
  return 0;
}
