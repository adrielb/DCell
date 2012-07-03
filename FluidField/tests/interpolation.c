#include "FluidField.h"

void InterfacialForceAdhesion(IrregularNode *n, void *context );

typedef struct _Context {
  PetscReal K;  // magnitude of surface tension
  PetscReal fa; // magnitude of adhesion
  PetscReal scale;
} Context;

void InterfacialForceAdhesion(IrregularNode *n, void *context )
{
//  const Context *c = (Context*)context;

//  n->F1 = c->scale*(-c->K * n->k);
  n->F1 = -n->k*100;
  n->F2 = 0;
}

int main(int argc, char **args) {
  PetscErrorCode ierr;
  ierr = DCellInit(); CHKERRQ(ierr);
  ierr = PetscInfoAllow(PETSC_TRUE,"info.log"); CHKERRQ(ierr);

  PetscReal d1 = 16, dx = 1./(d1-1);
  iCoor size = {d1,d1,0};

  FluidField fluid;
  ierr = FluidFieldCreate(PETSC_COMM_WORLD, &fluid);  CHKERRQ(ierr);
  ierr = FluidFieldSetDims(fluid,size); CHKERRQ(ierr);
  ierr = FluidFieldSetDx(fluid,dx); CHKERRQ(ierr);
  ierr = FluidFieldSetup(fluid); CHKERRQ(ierr);

  ierr = KSPSetTolerances(fluid->ksp,1e-12,1e-10,PETSC_DEFAULT,PETSC_DEFAULT); CHKERRQ(ierr);
  ierr = KSPMonitorSet(fluid->ksp,KSPMonitorDefault,0,0); CHKERRQ(ierr);

  IIM iim;
  ierr = IIMCreate( !fluid->is3D, 64, fluid->dh, &iim); CHKERRQ(ierr);
  ierr = IIMSetViscosity(iim, fluid->mu); CHKERRQ(ierr);

  ierr = IIMSetForceComponents(iim,InterfacialForceAdhesion ); CHKERRQ(ierr);

  LevelSet ls;
//  ierr = LevelSetInitializeToCircle(fluid->dh, (Coor) {0.5,0.5,0}, 0.2, &ls); CHKERRQ(ierr);
  ierr = LevelSetInitializeToStar2D(fluid->dh, (Coor) {0.5,0.5,0}, 0.3, 0.1, 3, &ls); CHKERRQ(ierr);

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

  int factor = 20;
  PetscReal di = 1.0 / factor;
  iCoor s = {size.x * factor, size.y * factor, 0};
  Grid g;
  ierr = GridCreate(fluid->dh, (iCoor){0,0,0}, s, 4, &g); CHKERRQ(ierr);

  iCoor p, q;
  int x, y;
  Coor  X;
  PetscReal ***vel, ***grid;
  ierr = DMDAVecGetArrayDOF(fluid->daV,fluid->vel,&vel); CHKERRQ(ierr);
  ierr = GridGet( g, &grid); CHKERRQ(ierr);
  ierr = GridGetBounds(g, &p, &q); CHKERRQ(ierr);
  for (y = p.y; y < q.y - 1.5*factor; y++ ) {
    for (x = p.x; x < q.x - 1.5*factor; x++ ) {
      Coor V = {0,0,0};
      X.x = x * di;
      X.y = y * di;
      ierr = InterpolateVelocity2D( 1, vel, X, &V ); CHKERRQ(ierr);
      grid[y][x][0] = V.x;
      grid[y][x][1] = V.y;
      ierr = IIMCorrectVelocity( iim, X, &V ); CHKERRQ(ierr);
      grid[y][x][2] = V.x;
      grid[y][x][3] = V.y;
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
