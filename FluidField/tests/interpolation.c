#include "FluidField.h"

void InterfacialForceAdhesion(IrregularNode *n, void *context );

typedef struct _Context {
  PetscReal K;  // magnitude of surface tension
  PetscReal fa; // magnitude of adhesion
  PetscReal scale;
} Context;

void InterfacialForceAdhesion(IrregularNode *n, void *context )
{
  const Context *c = (Context*)context;

//  n->f1 = c->scale*(-c->K * n->k);
  n->f1 = c->scale;
  n->f2 = c->scale;
}

int main(int argc, char **args) {
  PetscErrorCode ierr;
  ierr = DCellInit(); CHKERRQ(ierr);
  ierr = PetscInfoAllow(PETSC_TRUE,"info.log"); CHKERRQ(ierr);

  PetscReal d1 = 32, dx = 1./(d1-1);
  iCoor size = {d1,d1,0};

  FluidField fluid;
  ierr = FluidFieldCreate(PETSC_COMM_WORLD, &fluid);  CHKERRQ(ierr);
  ierr = FluidFieldSetDims(fluid,size); CHKERRQ(ierr);
  ierr = FluidFieldSetDx(fluid,dx); CHKERRQ(ierr);
  ierr = FluidFieldSetup(fluid); CHKERRQ(ierr);

  IIM iim;
  ierr = IIMCreate( !fluid->is3D, 64, fluid->dh, &iim); CHKERRQ(ierr);
  ierr = IIMSetViscosity(iim, fluid->mu); CHKERRQ(ierr);

  ierr = IIMSetForceComponents(iim,InterfacialForceAdhesion ); CHKERRQ(ierr);

  LevelSet ls;
  ierr = LevelSetInitializeToCircle(fluid->dh, (Coor) {0.5,0.5,0}, 0.4, &ls); CHKERRQ(ierr);

  Context context;
  context.K = 10;
  context.fa = 100;
  context.scale = 1e-3;

  ierr = IIMSetForceContext(iim, &context); CHKERRQ(ierr);
  ierr = GAPutVec(fluid->rhs,fluid->ga); CHKERRQ(ierr);
  ierr = IIMUpdateRHS(iim, ls, fluid->ga); CHKERRQ(ierr);
  ierr = GAGetVec(fluid->ga, fluid->rhs); CHKERRQ(ierr);
  ierr = FluidFieldSolve( fluid ); CHKERRQ(ierr);
  ierr = FluidFieldWrite( fluid, 0); CHKERRQ(ierr);
  ierr = LevelSetWriteIrregularNodeList(ls, 0); CHKERRQ(ierr);

  PetscReal di = 0.1;
  iCoor s = {fluid->lens.x / di, fluid->lens.y / di, 0};
  Grid g;
  ierr = GridCreate(fluid->dh, (iCoor){0,0,0}, s, 1, &g); CHKERRQ(ierr);

  int x, y;
  Coor  X;
  PetscReal ***vel, **grid;
  ierr = DMDAVecGetArrayDOF(fluid->daV,fluid->vel,&vel); CHKERRQ(ierr);
  ierr = GridGet( g, &grid); CHKERRQ(ierr);
  for (y = 0; y < s.y; y++ ) {
    for (x = 0; x < s.x; x++ ) {
      Coor V = {0,0,0};
      PetscReal *v = &V.x;
      X.x = x * di;
      X.y = y * di;
      ierr = InterpolateVelocity2D( 0, vel, X, &V ); CHKERRQ(ierr);
//      ierr = IIMCorrectVelocity( iim, X, &V ); CHKERRQ(ierr);
//      grid[y][x] = v[0];
      grid[y][x] = v[0];
    }
  }
  ierr = DMDAVecRestoreArrayDOF(fluid->daV,fluid->vel,&vel); CHKERRQ(ierr);

  ierr = GridWrite(g, 0); CHKERRQ(ierr);

  ierr = LevelSetDestroy(ls); CHKERRQ(ierr);
  ierr = IIMDestroy(iim); CHKERRQ(ierr);
  ierr = FluidFieldDestroy(fluid); CHKERRQ(ierr);
  ierr = DCellFinalize(); CHKERRQ(ierr);
  return 0;
}
