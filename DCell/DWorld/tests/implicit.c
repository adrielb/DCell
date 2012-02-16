#include "DWorld.h"

typedef struct _MyCell {
  struct _DCell dcell;
  PetscReal K;  // magnitude of surface tension
} *MyCell;

PetscErrorCode MyCellWrite( DCell dcell, int ti )
{
  PetscErrorCode ierr;
  PetscFunctionBegin;
  ierr = DCellWrite(dcell, ti); CHKERRQ(ierr);
  ierr = ParticleLSWriteParticles(dcell->lsPlasmaMembrane->pls, ti); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

inline PetscReal ClipCurvature( PetscReal k )
{
  const PetscReal clip = 0.1;

  k = k >  clip ?  clip : k;
  k = k < -clip ? -clip : k;

  return k;
}

void InterfacialForceSurfaceTension(IrregularNode *n, void *context )
{
  const MyCell c = (MyCell)context;

  n->f1 = -c->K * ClipCurvature(n->k);
  n->f2 = 0;
}

#undef __FUNCT__
#define __FUNCT__ "UpdateFluidFieldRHS"
PetscErrorCode UpdateFluidFieldRHS( DCell dcell, IIM iim, int ga, PetscReal t ) {
  PetscErrorCode ierr;
  PetscFunctionBegin;
  ierr = IIMSetForceComponents(iim, InterfacialForceSurfaceTension ); CHKERRQ(ierr);
  ierr = IIMSetForceContext(iim, dcell); CHKERRQ(ierr);
  ierr = IIMUpdateRHS(iim, dcell->lsPlasmaMembrane, ga); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MyCellUpdateFluidFieldImplicitRHS"
PetscErrorCode MyCellUpdateFluidFieldImplicitRHS( DCell dcell, IIM iim, int ga, PetscReal t ) {
  PetscErrorCode ierr;
  PetscFunctionBegin;
  ierr = IIMSetForceComponents(iim, InterfacialForceSurfaceTension ); CHKERRQ(ierr);
  ierr = IIMSetForceContext(iim, dcell); CHKERRQ(ierr);
  LevelSet ls = dcell->lsPlasmaMembrane->psi;
  ierr = LevelSetUpdateIrregularNodeList(ls); CHKERRQ(ierr);
  ierr = IIMUpdateRHS(iim, ls, ga); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MyCellCreate"
PetscErrorCode MyCellCreate( LevelSet ls, MyCell *mycell )
{
  MyCell cell;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscNew( struct _MyCell, &cell); CHKERRQ(ierr);
  ierr = DCellSetup( ls, (DCell)cell ); CHKERRQ(ierr);
  cell->dcell.Write = MyCellWrite;
  cell->dcell.UpdateFluidFieldRHS = MyCellUpdateFluidFieldImplicitRHS;
  cell->K = 20;
  *mycell = cell;
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc, char **args) {
  PetscErrorCode ierr;
  ierr = DCellInit(); CHKERRQ(ierr);

  PetscReal dx = 0.1;
  Coor len = {24,24,0};
  iCoor size = {len.x/dx,len.y/dx,0};
  printf("MX = %d;\n", size.x);
  printf("MY = %d;\n", size.y);

  FluidField fluid;
  ierr = FluidFieldCreate(PETSC_COMM_WORLD, &fluid);  CHKERRQ(ierr);
  ierr = FluidFieldSetDims(fluid,size); CHKERRQ(ierr);
  ierr = FluidFieldSetDx(fluid,dx); CHKERRQ(ierr);
  ierr = FluidFieldSetup(fluid); CHKERRQ(ierr);

  DWorld world;
  ierr = DWorldCreate(fluid, &world); CHKERRQ(ierr);

  Coor center = (Coor){12,12,0};
  PetscReal radius = 9;
  LevelSet ls;
  ierr = LevelSetInitializeToStar2D(fluid->dh,center,radius,0.5*radius, 5,&ls); CHKERRQ(ierr);
  ls->Advect = LevelSetAdvectImplicit;

  MyCell cell;
  ierr = MyCellCreate( ls, &cell ); CHKERRQ(ierr);

  ierr = MyCellWrite((DCell)cell, 0); CHKERRQ(ierr);

  ierr = DWorldAddDCell( world, cell ); CHKERRQ(ierr);
  world->timax = 100;
  world->dtmax = 0.0001;
  world->CFL = 0.5;
  world->tend = 3;

  ierr = DWorldSimulate_BFGS(world); CHKERRQ(ierr);
  ierr = DWorldDestroy(world); CHKERRQ(ierr);
  ierr = DCellFinalize(); CHKERRQ(ierr);
  return 0;
}
