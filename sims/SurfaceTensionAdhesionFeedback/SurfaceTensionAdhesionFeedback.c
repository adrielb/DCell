#include "DWorld.h"

PetscReal DiracDelta( PetscReal **phi, Coor dh, const PetscReal x, const PetscReal y );
void InterfacialForceAdhesion(IrregularNode *n, void *context );

typedef struct _MyCell {
  struct _DCell dcell;
  PetscReal K;  // magnitude of surface tension
  PetscReal Ka; // scale factor of contact area on surface tension
  PetscReal fa; // magnitude of adhesion
  PetscReal scale;
  PetscReal contactThres;
  PetscReal contactArea;
  Coor dh;
  PetscViewer contactareafile;
} *MyCell;

#undef __FUNCT__
#define __FUNCT__ "MyCellCreate"
PetscErrorCode MyCellCreate( LevelSet ls, MyCell *mycell )
{
  MyCell cell;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscNew( struct _MyCell, &cell); CHKERRQ(ierr);
  ierr = DCellSetup( ls, (DCell)cell ); CHKERRQ(ierr);

  *mycell = cell;
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "CalcContactArea"
PetscErrorCode CalcContactArea( MyCell cell )
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
  ierr = PetscViewerASCIIPrintf(cell->contactareafile,"%e\n",cell->contactArea); CHKERRQ(ierr);
  return 0;
}

void InterfacialForceAdhesionFeedback(IrregularNode *n, void *context )
{
  const MyCell c = (MyCell)context;

  if (n->X.y < c->contactThres) {
    n->fa1 = c->fa/c->dh.x; //direction of adhesion normal to membrane
    n->fa2 = 0;
  }

  const PetscReal Fk = (c->K + c->contactArea * c->fa*c->fa * c->Ka )* n->k;
  n->f1 = c->scale*(n->fa1 - Fk);
  n->f2 = c->scale*(n->fa2);
}

void ECMFunction( MyCell cell, PetscReal t ) {
  if( t > 0.175329 ) {
    cell->fa = 10;
    return;
  }
  if( t > 0.00230684 ) {
    cell->fa = 200;
    return;
  }
}

#undef __FUNCT__
#define __FUNCT__ "UpdateFluidFieldRHS"
PetscErrorCode UpdateFluidFieldRHS( DCell dcell, IIM iim, int ga, PetscReal t ) {
  MyCell mycell = (MyCell) dcell;
  PetscErrorCode ierr;
  PetscFunctionBegin
//  ECMFunction(mycell, t);
  CalcContactArea(mycell);
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
#define __FUNCT__ "main"
int main(int argc, char **args) {
  PetscErrorCode ierr;
  ierr = DCellInit(); CHKERRQ(ierr);

  int d1 = 64;
  PetscReal dx = 1./d1;
  FluidField fluid;
  ierr = FluidFieldCreate(PETSC_COMM_WORLD, &fluid);  CHKERRQ(ierr);
  ierr = FluidFieldSetDims(fluid,(iCoor){2*d1,d1}); CHKERRQ(ierr);
  ierr = FluidFieldSetDx(fluid,dx); CHKERRQ(ierr);
  ierr = FluidFieldSetup(fluid); CHKERRQ(ierr);

  DWorld world;
  ierr = DWorldCreate(fluid, &world); CHKERRQ(ierr);

  Coor center = (Coor){1,.4,0};
  PetscReal radius = 0.4-dx/2.;
  LevelSet ls;
  ierr = LevelSetInitializeToCircle(fluid->dh,center,radius,&ls); CHKERRQ(ierr);
  MyCell cell;
  ierr = MyCellCreate( ls, &cell ); CHKERRQ(ierr);
  cell->dh = fluid->dh;
  cell->Ka= 0.002;
  cell->K = 10;
  cell->fa = 320;
  cell->scale = 1e-9;
  cell->contactThres = 1.5;
  cell->dcell.UpdateFluidFieldRHS = UpdateFluidFieldRHS;
  cell->dcell.Destroy = MyCellDestroy;
    char filename[PETSC_MAX_PATH_LEN];
    char wd[PETSC_MAX_PATH_LEN];
    ierr = DCellGetWorkingDirectory(wd); CHKERRQ(ierr);
    ierr = PetscSNPrintf(filename,PETSC_MAX_PATH_LEN,"%s/contactarea.dat",wd); CHKERRQ(ierr);
    ierr = PetscViewerASCIIOpen(PETSC_COMM_SELF,filename,&cell->contactareafile); CHKERRQ(ierr);
  ierr = DWorldAddDCell( world, cell ); CHKERRQ(ierr);
  world->timax = 200;
  world->dtmax = 5e-1;
  world->CFL = 0.9;
  world->tend = 1e3;

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

  if (n->X.y < c->contactThres && n->ny < -0.707) {
    n->fa1 = c->fa; //direction of adhesion normal to membrane
    n->fa2 = 0;
  }

  n->f1 = c->scale*(n->fa1 - c->K * n->k);
  n->f2 = c->scale*(n->fa2);
}
