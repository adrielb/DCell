#include "DWorld.h"

PetscReal DiracDelta( PetscReal **phi, Coor dh, const PetscReal x, const PetscReal y );

typedef struct _MyCell {
  struct _DCell dcell;
  PetscReal K;  // magnitude of surface tension
  PetscReal Ka; // scale factor of contact area on surface tension
  PetscReal scale;
  PetscReal contactThres;
  PetscReal contactArea;
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
#define __FUNCT__ "MyCellCreate"
PetscErrorCode MyCellCreate( LevelSet ls, MyCell *mycell )
{
  MyCell cell;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscNew( struct _MyCell, &cell); CHKERRQ(ierr);
  ierr = DCellSetup( ls, (DCell)cell ); CHKERRQ(ierr);
  cell->dcell.Write = MyCellWrite;
  *mycell = cell;
  PetscFunctionReturn(0);
}

PetscReal AdhesionForce( PetscReal x ) {
  const PetscReal b = 100;
  const PetscReal a = 0;
  return b + a * x;
}

void InterfacialForceAdhesion(IrregularNode *n, void *context )
{
  const MyCell c = (MyCell)context;

  if (n->X.y < c->contactThres && n->ny < 0) {
    n->fa1 = AdhesionForce( n->X.x * c->dh.x ) / c->dh.x; //direction of adhesion normal to membrane
    n->fa2 = 0;
  }
  n->f1 = c->scale*(n->fa1 - c->K * n->k);
  n->f2 = c->scale*(n->fa2);
}

PetscErrorCode UpdateFluidFieldRHS( DCell dcell, IIM iim, int ga, PetscReal t ) {
  PetscErrorCode ierr;
  PetscFunctionBegin
  ierr = IIMSetForceComponents(iim, InterfacialForceAdhesion ); CHKERRQ(ierr);
  ierr = IIMSetForceContext(iim, dcell); CHKERRQ(ierr);
  ierr = IIMUpdateRHS(iim, dcell->lsPlasmaMembrane, ga); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

int main(int argc, char **args) {
  PetscErrorCode ierr;
  ierr = DCellInit(); CHKERRQ(ierr);

  PetscReal dx = 1;
  Coor len = {50,25,0};
  Coor dh  = {dx,dx,0};
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

  Coor center = (Coor){10,10,0};
  PetscReal radius = 10-dx;
  LevelSet ls;
  ierr = LevelSetInitializeToCircle(fluid->dh,center,radius,&ls); CHKERRQ(ierr);
  ierr = LevelSetInitializeParticles(ls); CHKERRQ(ierr);
  MyCell cell;
  ierr = MyCellCreate( ls, &cell ); CHKERRQ(ierr);
  cell->dh = fluid->dh;
  cell->Ka= 0.002;
  cell->K = 200;
  cell->scale = 1e-6;
  cell->contactThres = 1.5;
  cell->dcell.UpdateFluidFieldRHS = UpdateFluidFieldRHS;
  ierr = DWorldAddDCell( world, cell ); CHKERRQ(ierr);
  world->timax = 400;
  world->dtmax = 12000;
  world->CFL = 0.9;
  world->tend = 1e30;

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

/*
 * M1 F11  Debug
 * M2 F8   Resume
 * M3 F7   Step return
 * M4 F6   Step over
 * M5 F5   Step into
 */
