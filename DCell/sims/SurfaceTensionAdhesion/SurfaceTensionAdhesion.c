#include "DWorld.h"

PetscReal DiracDelta( PetscReal **phi, Coor dh, const PetscReal x, const PetscReal y );
void InterfacialForceAdhesion(IrregularNode *n, void *context );

typedef struct _MyCell {
  struct _DCell dcell;
  PetscReal k0; // original curvature of undeformed sphere
  PetscReal K;  // magnitude of surface tension
  PetscReal fa; // magnitude of adhesion
  PetscReal scale;
  PetscReal contactThres;
  PetscReal contactArea;
  Coor dh;
  PetscViewer contactareafile;
} *MyCell;

#undef __FUNCT__
#define __FUNCT__ "CalcContactArea"
PetscErrorCode CalcContactArea( MyCell cell, PetscReal t )
{
  LevelSet ls = cell->dcell.lsPlasmaMembrane;
  PetscReal **phi;
  iCoor p,q;
  PetscReal x,y;
  PetscReal dh = 0.1;
  PetscReal dA = dh*dh*ls->phi->d.x*ls->phi->d.y;
  PetscErrorCode ierr;

  ierr = GridGetBounds(ls->phi,&p,&q); CHKERRQ(ierr);
  ierr = GridGet(ls->phi,&phi); CHKERRQ(ierr);
  cell->contactArea = 0;
  for (y = p.y; y < q.y && y < cell->contactThres; y += dh ) {
    for (x = p.x; x < q.x; x += dh ) {
      cell->contactArea += DiracDelta( phi, ls->phi->d, x, y )*dA;
    }
  }
  ierr = PetscViewerASCIIPrintf(cell->contactareafile,"%e %e\n", t, cell->contactArea); CHKERRQ(ierr);
  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "UpdateFluidFieldRHS"
PetscErrorCode UpdateFluidFieldRHS( DCell dcell, IIM iim, int ga, PetscReal t ) {
  MyCell mycell = (MyCell) dcell;
  PetscErrorCode ierr;
  PetscFunctionBegin
//  ECMFunction(mycell, t);
  CalcContactArea(mycell, t);
  ierr = IIMSetForceComponents(iim, InterfacialForceAdhesion ); CHKERRQ(ierr);
  ierr = IIMSetForceContext(iim, dcell); CHKERRQ(ierr);
  ierr = IIMUpdateRHS(iim, dcell->lsPlasmaMembrane, ga); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MyCellDestroy"
PetscErrorCode MyCellDestroy( DCell cell ) {
  MyCell mycell = (MyCell) cell;
  PetscErrorCode ierr;
  PetscFunctionBegin
  ierr = PetscViewerDestroy(mycell->contactareafile); CHKERRQ(ierr);
  ierr = DCellDestroy(cell); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MyCellWrite"
PetscErrorCode MyCellWrite( DCell dcell, int ti ) {
  PetscErrorCode ierr;
  PetscFunctionBegin
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
  cell->dcell.UpdateFluidFieldRHS = UpdateFluidFieldRHS;
  cell->dcell.Destroy = MyCellDestroy;
  cell->dcell.Write = MyCellWrite;
  *mycell = cell;
  ierr = PetscInfo(0,"Created MyCell\n"); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc, char **args) {
  PetscErrorCode ierr;
  ierr = DCellInit(); CHKERRQ(ierr);

  PetscReal dx = 0.5;
  Coor len = {50,25,0};
  iCoor size = {len.x/dx,len.y/dx,0};
  printf("Domain Size: %d x %d\n", size.x, size.y);

  FluidField fluid;
  ierr = FluidFieldCreate(PETSC_COMM_WORLD, &fluid);  CHKERRQ(ierr);
  ierr = FluidFieldSetDims(fluid,size); CHKERRQ(ierr);
  ierr = FluidFieldSetDx(fluid,dx); CHKERRQ(ierr);
  fluid->mu = 10;
  ierr = FluidFieldSetup(fluid); CHKERRQ(ierr);

  DWorld world;
  ierr = DWorldCreate(fluid, &world); CHKERRQ(ierr);

  PetscReal radius = 10;
  Coor center = (Coor){len.x/2.,radius+dx/2,0};
  LevelSet ls;
  ierr = LevelSetInitializeToCircle(fluid->dh,center,radius,&ls); CHKERRQ(ierr);
  ierr = LevelSetInitializeParticles(ls); CHKERRQ(ierr);

  MyCell cell;
  ierr = MyCellCreate( ls, &cell ); CHKERRQ(ierr);
  cell->k0 = 1 / radius;
  cell->dh = fluid->dh;
  cell->K = 50;
  cell->fa = 50;
  cell->scale = 1;
  cell->contactThres = 3.1;

    ierr = PetscOptionsGetReal(0,"-fa",&cell->fa,0); CHKERRQ(ierr);
    ierr = PetscOptionsGetReal(0,"-K",&cell->K,0); CHKERRQ(ierr);
    char filename[PETSC_MAX_PATH_LEN];
    char wd[PETSC_MAX_PATH_LEN];
    ierr = DCellGetWorkingDirectory(wd); CHKERRQ(ierr);
    ierr = PetscSNPrintf(filename,PETSC_MAX_PATH_LEN,"%s/contactarea.%.0f.%.0f.dat", wd, cell->fa, cell->K); CHKERRQ(ierr);
    ierr = PetscViewerASCIIOpen(PETSC_COMM_SELF,filename,&cell->contactareafile); CHKERRQ(ierr);

  ierr = DWorldAddDCell( world, cell ); CHKERRQ(ierr);
  world->timax = 2000;
  world->dtmax = 0.4;
  world->CFL = 0.1;
  world->tend = 8e2;
  world->dtframe = world->tend / 800.;
//  world->writeInterval = 10;

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

void InterfacialForceAdhesion(IrregularNode *n, void *context )
{
  const MyCell c = (MyCell)context;
  PetscReal K = c->K, ny = n->ny, k0 = c->k0;

  /*
  double clip = 0.1, k = n->k;
  k = k >  clip ?  clip : k;
  k = k < -clip ? -clip : k;
  */

  Coor b = {0,-1,0}; // Normal direction of contact surface.
  const PetscReal dot = n->nx * b.x + n->ny * b.y;
  if ( n->X.y < c->contactThres ) { // && n->ny < -0.707
//    n->fa1 = c->fa * (c->contactThres - n->X.y) / c->contactThres; //direction of adhesion normal to membrane
//    ny = ny > 0 ? 0 : ny;

    /* normal to flat surface
    n->fa1 = c->fa * ( n->nx * b.x + n->ny * b.y );
    n->fa2 = c->fa * ( n->sx * b.x + n->sy * b.y );
    K = 2 * c->K;
    */

    k0 = 1;
    if( dot > 0 ) {
      // normal to cell surface
      n->fa1 = c->fa;
      n->fa2 = 0;
    }
  }

  n->f1 = c->scale*(n->fa1 - K * ( n->k_nn - k0 ) );
  n->f2 = c->scale*(n->fa2);
}
