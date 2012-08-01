#include "DWorld.h"

PetscReal DiracDelta( PetscReal **phi, Coor dh, const PetscReal x, const PetscReal y );

typedef struct _MyCell {
  struct _DCell dcell;
  PetscReal k0;
  PetscReal K;  // magnitude of surface tension
  Coor dh;
} *MyCell;

PetscErrorCode MyCellWrite( DCell dcell, int ti )
{
  PetscErrorCode ierr;
  PetscFunctionBegin;
  ierr = DCellWrite(dcell, ti); CHKERRQ(ierr);
  ierr = ParticleLSWriteParticles(dcell->lsPlasmaMembrane->pls, ti); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "MyCellUpdateFluidFieldRHS"
PetscErrorCode MyCellUpdateFluidFieldRHS( DCell dcell, IIM iim, int ga, PetscReal t ) {
  PetscErrorCode ierr;
  PetscFunctionBegin;
  ierr = IIMSetForceComponents(iim, InterfacialForceSurfaceTension ); CHKERRQ(ierr);
  ierr = IIMSetForceContext(iim, dcell); CHKERRQ(ierr);
  ierr = IIMUpdateRHS(iim, dcell->lsPlasmaMembrane, ga); CHKERRQ(ierr);
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
  cell->dcell.UpdateFluidFieldRHS = MyCellUpdateFluidFieldRHS;
  *mycell = cell;
  PetscFunctionReturn(0);
}

void InterfacialForceSurfaceTension(IrregularNode *n, void *context )
{
  const MyCell c = (MyCell)context;

  n->F1 = -10;
  n->F1 = -c->K * n->k;
//  n->F1 = -(k - 10*c->k0);
  n->f2 = 0;
}

inline PetscReal ClipCurvature( PetscReal k )
{
  const PetscReal clip = 0.2;

  k = k >  clip ?  clip : k;
  k = k < -clip ? -clip : k;

  return k;
}

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc, char **args) {
  PetscErrorCode ierr;
  ierr = DCellInit(); CHKERRQ(ierr);

  FluidField fluid;
  ierr = FluidFieldCreate(PETSC_COMM_WORLD, &fluid);  CHKERRQ(ierr);
  ierr = FluidFieldSetup(fluid); CHKERRQ(ierr);

  DWorld world;
  ierr = DWorldCreate(fluid, &world); CHKERRQ(ierr);

  Coor center = (Coor){12,12,0};
  PetscReal radius = 9;
  LevelSet ls;
  ierr = LevelSetInitializeToStar2D(fluid->dh,center,radius,0.5*radius, 5,&ls); CHKERRQ(ierr);
//  ierr = LevelSetInitializeToCircle(fluid->dh, center, radius, &ls); CHKERRQ(ierr);
  ierr = LevelSetInitializeParticles(ls); CHKERRQ(ierr);
  MyCell cell;
  ierr = MyCellCreate( ls, &cell ); CHKERRQ(ierr);

  ierr = MyCellWrite((DCell)cell, 0); CHKERRQ(ierr);
  cell->k0 = 1 / radius;
  cell->dh = fluid->dh;
  cell->K = 20;

  ierr = DWorldAddDCell( world, cell ); CHKERRQ(ierr);
  ierr = DWorldSetFromOptions(world); CHKERRQ(ierr);
  world->Simulate = DWorldSimulate_Euler;
  ierr = DWorldSimulate(world); CHKERRQ(ierr);
  ierr = DWorldDestroy(world); CHKERRQ(ierr);
  ierr = DCellFinalize(); CHKERRQ(ierr);
  return 0;
}

inline PetscReal DiracDelta( PetscReal **phi, Coor dh, const PetscReal x, const PetscReal y )
{
  const PetscReal eps = 1.5;
  const PetscReal phival = Bilinear2D(GridFunction2D_Identity,phi,dh,x,y);
  if( phival < -eps || phival > eps ) return 0;
  const PetscReal px = Bilinear2D(GridFunction2D_DerivX,phi,dh,x,y);
  const PetscReal py = Bilinear2D(GridFunction2D_DerivY,phi,dh,x,y);
  const PetscReal gradmag = PetscSqrtScalar( px*px + py*py );
  return gradmag * ( 1 + cos( PETSC_PI * phival / eps) ) / (2*eps);
}
