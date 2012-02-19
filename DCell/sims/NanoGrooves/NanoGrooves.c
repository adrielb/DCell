#include "DWorld.h"

static LevelSet lsGrooves;

typedef struct _MyCell {
  struct _DCell dcell;
  PetscReal Fk;  // magnitude of surface tension
  PetscReal Fk0;  // magnitude of basal surface tension
  PetscReal Fa;  // magnitude of adhesion
  PetscReal kclip; // limit of curvature force
  PetscReal ecm; // ecm concentration
  PetscReal scale;
  PetscReal contactThres;
  PetscReal contactArea;
  Coor dh;
  PetscViewer contactareafile;
} *MyCell;

PetscErrorCode CalcContactArea( MyCell cell, PetscReal t );
PetscErrorCode MyCellSetFromOptions( MyCell cell );
PetscErrorCode MyCellCreate( LevelSet ls, MyCell *mycell );
PetscErrorCode MyCellUpdateFluidFieldRHS( DCell dcell, IIM iim, int ga, PetscReal t );
PetscErrorCode MyCellUpdateFluidFieldImplicitRHS( DCell dcell, IIM iim, int ga, PetscReal t );
PetscErrorCode MyCellWrite( DCell dcell, int ti );

void InterfacialForceAdhesion(IrregularNode *n, void *context )
{
  const MyCell c = (MyCell)context;
  const Coor dh = lsGrooves->phi->d;
  PetscReal dist;
  PetscReal ecmTot = c->ecm * c->contactArea;

  GridInterpolate( lsGrooves->phi, n->X, &dist );

  if( dist*dh.x > -c->contactThres ) {
    n->fa1 = c->Fa;
    n->fa2 = 0;
  }

  PetscReal k = n->k;
  const PetscReal clip = c->kclip;
  k = k >  clip ?  clip : k;
  k = k < -clip ? -clip : k;

  n->f1 = c->scale * ( n->fa1*c->ecm - (c->Fk * ecmTot + c->Fk0 ) * k);
  n->f2 = c->scale * ( n->fa2 );
}

int main(int argc, char **args) {
  PetscErrorCode ierr;
  ierr = DCellInit(); CHKERRQ(ierr);

  PetscReal dx = 0.1;
  Coor len = {100,25,0};
  Coor dh = {dx,dx,0};
  iCoor size = {len.x/dx,len.y/dx,0};
  printf("MX = %d;\n", size.x);
  printf("MY = %d;\n", size.y);

  iCoor pos = {0,0,0};
  ierr = LevelSetCreate(dh,pos,size,&lsGrooves); CHKERRQ(ierr);
  ierr = VecSet(lsGrooves->phi->v,-1); CHKERRQ(ierr);
  int borderwidth = 3;
  ierr = GridDrawBorder( lsGrooves->phi, borderwidth, 1 ); CHKERRQ(ierr);
  //Draw posts
  PetscReal c,
            w = 1.0,
            h = 0.5+borderwidth*dx,
            s = 2.0;

  ierr = PetscOptionsGetReal(0,"-groove_width",&w,0); CHKERRQ(ierr);

  for ( c = 0; c < len.x; c+=s ) {
    Coor lo = {c,  0,0};
    Coor hi = {c+w,h,0};
    ierr = GridFillRectangle(lsGrooves->phi, lo, hi, 1); CHKERRQ(ierr);
  }
  ierr = GridSetName(lsGrooves->phi,"grooves"); CHKERRQ(ierr);
  ierr = GridWrite(lsGrooves->phi,0); CHKERRQ(ierr);
  ierr = LevelSetInitializeFromImage(lsGrooves); CHKERRQ(ierr);
  ierr = GridWrite(lsGrooves->phi,1); CHKERRQ(ierr);
  ierr = LevelSetWriteIrregularNodeList(lsGrooves,0); CHKERRQ(ierr);

  FluidField fluid;
  ierr = FluidFieldCreate(PETSC_COMM_WORLD, &fluid);  CHKERRQ(ierr);
  ierr = FluidFieldSetDims(fluid,size); CHKERRQ(ierr);
  ierr = FluidFieldSetDx(fluid,dx); CHKERRQ(ierr);
  ierr = FluidFieldSetMask(fluid, lsGrooves->phi ); CHKERRQ(ierr);
  ierr = FluidFieldSetup(fluid); CHKERRQ(ierr);

  DWorld world;
  ierr = DWorldCreate(fluid, &world); CHKERRQ(ierr);
  ierr = DWorldSetPrintStep(world,PETSC_TRUE); CHKERRQ(ierr);

  PetscReal radius = 10;
  Coor center = (Coor){len.x/2,radius+h,0};
  LevelSet ls;
  ierr = LevelSetInitializeToCircle(fluid->dh,center,radius,&ls); CHKERRQ(ierr);
//  ierr = LevelSetInitializeParticles(ls); CHKERRQ(ierr);
  ls->Advect = LevelSetAdvectSL;
  ierr = GridSetName(ls->phi,"mycell"); CHKERRQ(ierr);

  MyCell cell;
  ierr = MyCellCreate( ls, &cell ); CHKERRQ(ierr);
  cell->dh = fluid->dh;
  cell->Fk = 3;
  cell->Fk = 1;
  cell->Fa = 3;
  cell->kclip = 1 / radius;
  cell->ecm = 1;
  cell->scale = 1e-3;
  cell->contactThres = 0.3;
  ierr = MyCellSetFromOptions(cell); CHKERRQ(ierr);

  char filename[PETSC_MAX_PATH_LEN];
  char wd[PETSC_MAX_PATH_LEN];
  ierr = DCellGetWorkingDirectory(wd); CHKERRQ(ierr);
  ierr = PetscSNPrintf(filename,PETSC_MAX_PATH_LEN,"%s/contactarea.dat", wd); CHKERRQ(ierr);
  ierr = PetscViewerASCIIOpen(PETSC_COMM_SELF,filename,&cell->contactareafile); CHKERRQ(ierr);

  ierr = DWorldAddDCell( world, cell ); CHKERRQ(ierr);
  world->timax = 10;
  world->dtmax = 1e9;
  world->CFL = 0.1;
  world->tend = 1e9;
  ierr = DWorldSetFromOptions(world); CHKERRQ(ierr);
  ierr = DWorldSimulate(world); CHKERRQ(ierr);

  ierr = LevelSetDestroy(lsGrooves); CHKERRQ(ierr);
  ierr = DWorldDestroy(world); CHKERRQ(ierr);
  ierr = DCellFinalize(); CHKERRQ(ierr);
  return 0;
}


#undef __FUNCT__
#define __FUNCT__ "MyCellDestroy"
PetscErrorCode MyCellDestroy( DCell cell ) {
  MyCell mycell = (MyCell) cell;
  PetscErrorCode ierr;
  PetscFunctionBegin;
  ierr = PetscViewerDestroy(&mycell->contactareafile); CHKERRQ(ierr);
  ierr = DCellDestroy(cell); CHKERRQ(ierr);
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
  cell->dcell.UpdateFluidFieldRHS = MyCellUpdateFluidFieldImplicitRHS;
  cell->dcell.Destroy = MyCellDestroy;
  cell->dcell.Write = MyCellWrite;
  *mycell = cell;
  ierr = PetscInfo(0,"Created MyCell\n"); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MyCellSetFromOptions"
PetscErrorCode MyCellSetFromOptions( MyCell cell )
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscOptionsGetReal(0,"-Fa",&cell->Fa,0); CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(0,"-Fk", &cell->Fk, 0); CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(0,"-Fk0", &cell->Fk0, 0); CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(0,"-kclip",&cell->kclip,0); CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(0,"-ecm",&cell->ecm,0); CHKERRQ(ierr);

  ierr = PetscPrintf(MPI_COMM_WORLD,"Cell Parameters\n"); CHKERRQ(ierr);
  ierr = PetscPrintf(MPI_COMM_WORLD,"Fa     = %f\n", cell->Fa); CHKERRQ(ierr);
  ierr = PetscPrintf(MPI_COMM_WORLD,"Fk     = %f\n", cell->Fk); CHKERRQ(ierr);
  ierr = PetscPrintf(MPI_COMM_WORLD,"Fk0    = %f\n", cell->Fk0); CHKERRQ(ierr);
  ierr = PetscPrintf(MPI_COMM_WORLD,"kclip  = %f\n", cell->kclip); CHKERRQ(ierr);
  ierr = PetscPrintf(MPI_COMM_WORLD,"ecm    = %f\n", cell->ecm); CHKERRQ(ierr);
  ierr = PetscPrintf(MPI_COMM_WORLD,"---------------\n", cell->ecm); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MyCellUpdateFluidFieldRHS"
PetscErrorCode MyCellUpdateFluidFieldRHS( DCell dcell, IIM iim, int ga, PetscReal t ) {
  PetscErrorCode ierr;
  PetscFunctionBegin;
  ierr = CalcContactArea( (MyCell)dcell, t ); CHKERRQ(ierr);
  ierr = IIMSetForceComponents(iim, InterfacialForceAdhesion ); CHKERRQ(ierr);
  ierr = IIMSetForceContext(iim, dcell); CHKERRQ(ierr);
  ierr = IIMUpdateRHS(iim, dcell->lsPlasmaMembrane, ga); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MyCellUpdateFluidFieldImplicitRHS"
PetscErrorCode MyCellUpdateFluidFieldImplicitRHS( DCell dcell, IIM iim, int ga, PetscReal t ) {
  PetscErrorCode ierr;
  PetscFunctionBegin;
  ierr = CalcContactArea( (MyCell)dcell, t ); CHKERRQ(ierr);
  ierr = IIMSetForceComponents(iim, InterfacialForceAdhesion ); CHKERRQ(ierr);
  ierr = IIMSetForceContext(iim, dcell); CHKERRQ(ierr);
  LevelSet ls = dcell->lsPlasmaMembrane->psi;
  ierr = LevelSetUpdateIrregularNodeList(ls); CHKERRQ(ierr);
  ierr = IIMUpdateRHS(iim, ls, ga); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MyCellWrite"
PetscErrorCode MyCellWrite( DCell dcell, int ti ) {
  PetscErrorCode ierr;
  PetscFunctionBegin;
  ierr = DCellWrite(dcell, ti); CHKERRQ(ierr);
  ierr = ParticleLSWriteParticles(dcell->lsPlasmaMembrane->pls, ti); CHKERRQ(ierr);
  ierr = LevelSetWriteIrregularNodeList(dcell->lsPlasmaMembrane->psi, ti); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "CalcContactArea"
PetscErrorCode CalcContactArea( MyCell cell, PetscReal t )
{
  LevelSet ls = cell->dcell.lsPlasmaMembrane;
  PetscReal **phi;
  iCoor p,q;
  Coor X;
  const int idh = 10; // samples per grid cell
  const PetscReal dh = 1. / idh;
  const PetscReal dA = dh*dh*ls->phi->d.x*ls->phi->d.y;
  const PetscReal dx = ls->phi->d.x;
  PetscReal dist;
  int m, i, j;
  const int len = ArrayLength(ls->band);
  const iCoor *band = ArrayGetData(ls->band);
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = GridGetBounds(ls->phi,&p,&q); CHKERRQ(ierr);
  ierr = GridGet(ls->phi,&phi); CHKERRQ(ierr);
  cell->contactArea = 0;
  for (m = 0; m < len; ++m) {
    for (j = 0; j < idh; j++) {
      for (i = 0; i < idh; i++) {
        X.x = band[m].x + i * dh;
        X.y = band[m].y + j * dh;
        ierr = GridInterpolate( lsGrooves->phi, X, &dist ); CHKERRQ(ierr);
        if( dist*dx > -cell->contactThres ) {
          cell->contactArea += LevelSetDiracDelta2D( phi, ls->phi->d, X )*dA;
        } // if in contact
      } //j
    } //i
  } //m in band
  ierr = PetscViewerASCIIPrintf(cell->contactareafile,"%e %e\n", t, cell->contactArea); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
