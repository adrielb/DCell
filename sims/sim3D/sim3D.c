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
  PetscReal adhesionRadius;
  Coor dh;
  PetscViewer contactareafile;
  SpatialIndex sidx;
} *MyCell;

PetscErrorCode MyCellCreate( LevelSet ls, MyCell *mycell );
PetscErrorCode MyCellDestroy( DCell cell );
PetscErrorCode MyCellSetFromOptions( MyCell cell );
PetscErrorCode MyCellUpdateFluidFieldRHS( DCell dcell, IIM iim, int ga, PetscReal t );
PetscErrorCode MyCellUpdateFluidFieldImplicitRHS( DCell dcell, IIM iim, int ga, PetscReal t );
PetscErrorCode MyCellWrite( DCell dcell, int ti );
PetscErrorCode MyCellWriteImplicit( DCell dcell, int ti );

void InterfacialForceAdhesion(IrregularNode *n, void *context )
{
  const MyCell c = (MyCell)context;
  const Coor dh = lsGrooves->phi->d;
  const PetscReal r = c->adhesionRadius;
  PetscReal dist;
  const int MAXLEN = 500;
  int len = 0;
  int i;
  IrregularNode *nei[MAXLEN];

  GridInterpolate( lsGrooves->phi, n->X, &dist );
  if( dist*dh.x > -c->contactThres ) {
    SpatialIndexQueryPoints( c->sidx, n->X, r / dh.x, MAXLEN, &len, (void**)nei);

    for (i = 0; i < len; ++i) {
      nei[i]->fa1 = c->Fa;
    }
  }

  PetscReal k = n->k;
  const PetscReal clip = c->kclip;
  k = k >  clip ?  clip : k;
  k = k < -clip ? -clip : k;

  n->F1 = c->scale * ( n->fa1 - c->Fk0 * k);
  n->F2 = c->scale * ( n->fa2 );
}

int main(int argc, char **args) {
  PetscErrorCode ierr;
  ierr = DCellInit(); CHKERRQ(ierr);

  FluidField fluid;
  ierr = FluidFieldCreate(PETSC_COMM_WORLD, &fluid);  CHKERRQ(ierr);

  Coor len = fluid->lens;
  Coor  dh = fluid->dh;
  iCoor size = fluid->dims;
  PetscReal dx = dh.x;

  // Make NanoSurface
  int borderwidth = 2;
  PetscReal c,
            b = borderwidth*dx,
            w = 1,
            h = 0.5,
            l = len.y,
            s = 2;
  ierr = PetscOptionsGetReal(0,"-groove_width",&w,0); CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(0,"-groove_height",&h,0); CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(0,"-groove_spacing",&s,0); CHKERRQ(ierr);
  iCoor pos = {0,0,0};
  ierr = LevelSetCreate(dh,pos,size,&lsGrooves); CHKERRQ(ierr);
  ierr = VecSet(lsGrooves->phi->v,-1); CHKERRQ(ierr);
  ierr = GridDrawBorder( lsGrooves->phi, borderwidth, 1 ); CHKERRQ(ierr);
  for ( c = dx-w/2; c < len.x; c+=s ) {
    Coor lo = {c,   0, 0};
    Coor hi = {c+w, l, h+b};
    ierr = GridFillBox(lsGrooves->phi, lo, hi, 1); CHKERRQ(ierr);
  }
  ierr = GridSetName(lsGrooves->phi,"grooves"); CHKERRQ(ierr);
  ierr = GridWrite(lsGrooves->phi,0); CHKERRQ(ierr);
  ierr = LevelSetInitializeFromImage(lsGrooves); CHKERRQ(ierr);
  ierr = GridWrite(lsGrooves->phi,1); CHKERRQ(ierr);
  ierr = LevelSetWriteIrregularNodeList(lsGrooves,0); CHKERRQ(ierr);


  ierr = FluidFieldSetMask(fluid, lsGrooves->phi ); CHKERRQ(ierr);
  ierr = FluidFieldSetup(fluid); CHKERRQ(ierr);

  DWorld world;
  ierr = DWorldCreate(fluid, &world); CHKERRQ(ierr);
  ierr = DWorldSetPrintStep(world,PETSC_TRUE); CHKERRQ(ierr);

  int rank;
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
  if( rank == 0 ) {
    PetscReal radius = 3;
    Coor center = (Coor){ len.x / 2., len.y / 2., radius + h + b - dx };
    LevelSet ls;
    ierr = LevelSetInitializeToSphere(fluid->dh,center,radius,&ls); CHKERRQ(ierr);
//  ierr = LevelSetInitializeToStar3D(fluid->dh, center, radius, 2.5, 5, &ls); CHKERRQ(ierr);
//  ierr = LevelSetInitializeParticles(ls); CHKERRQ(ierr);
    ierr = GridSetName(ls->phi,"mycell"); CHKERRQ(ierr);

    MyCell cell;
    ierr = MyCellCreate( ls, &cell ); CHKERRQ(ierr);
    cell->dh = fluid->dh;
    cell->kclip = 1 / radius;
    ierr = MyCellSetFromOptions(cell); CHKERRQ(ierr);
    ierr = DWorldAddDCell( world, cell ); CHKERRQ(ierr);

    ierr = VecShift(lsGrooves->phi->v,cell->contactThres/dx); CHKERRQ(ierr);
    ierr = LevelSetUpdateIrregularNodeList(lsGrooves); CHKERRQ(ierr);
    ierr = LevelSetWriteIrregularNodeList(lsGrooves,1); CHKERRQ(ierr);
    ierr = VecShift(lsGrooves->phi->v,-cell->contactThres/dx); CHKERRQ(ierr);
  }

  world->timax = 10;
  world->dtmax = 1e9;
  world->CFL = 0.1;
  world->tend = 1e9;
  ierr = DWorldSetFromOptions(world); CHKERRQ(ierr);
  world->Simulate = DWorldSimulate_Euler;
  ierr = DWorldSimulate(world); CHKERRQ(ierr);

  ierr = DWorldDestroy(world); CHKERRQ(ierr);
  ierr = DCellFinalize(); CHKERRQ(ierr);
  return 0;
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
  const PetscBool implicit = PETSC_FALSE;
  if( implicit ) {
    cell->dcell.UpdateFluidFieldRHS = MyCellUpdateFluidFieldImplicitRHS;
    cell->dcell.Write = MyCellWriteImplicit;
  } else {
    cell->dcell.UpdateFluidFieldRHS = MyCellUpdateFluidFieldRHS;
    cell->dcell.Write = MyCellWrite;
  }

  cell->dcell.Destroy = MyCellDestroy;

  cell->Fk = 1;
  cell->Fk0 = 1;
  cell->Fa = 1;
  cell->ecm = 1;
  cell->scale = 1e-3;
  cell->contactThres = 0.4;
  cell->adhesionRadius = 0.2;

  char filename[PETSC_MAX_PATH_LEN];
  char wd[PETSC_MAX_PATH_LEN];
  ierr = DCellGetWorkingDirectory(wd); CHKERRQ(ierr);
  ierr = PetscSNPrintf(filename,PETSC_MAX_PATH_LEN,"%s/contactarea.dat", wd); CHKERRQ(ierr);
  ierr = PetscViewerASCIIOpen(PETSC_COMM_SELF,filename,&cell->contactareafile); CHKERRQ(ierr);

  *mycell = cell;
  ierr = PetscInfo(0,"Created MyCell\n"); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MyCellDestroy"
PetscErrorCode MyCellDestroy( DCell cell )
{
  MyCell mycell = (MyCell) cell;
  PetscErrorCode ierr;
  PetscFunctionBegin;
  ierr = PetscViewerDestroy(&mycell->contactareafile); CHKERRQ(ierr);
  ierr = DCellDestroy(cell); CHKERRQ(ierr);
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
  ierr = PetscOptionsGetReal(0,"-contactThres",&cell->contactThres,0); CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(0,"-adhesionRadius",&cell->adhesionRadius,0); CHKERRQ(ierr);

  ierr = PetscPrintf(MPI_COMM_WORLD,"Cell Parameters\n"); CHKERRQ(ierr);
  ierr = PetscPrintf(MPI_COMM_WORLD,"Fa     = %f\n", cell->Fa); CHKERRQ(ierr);
  ierr = PetscPrintf(MPI_COMM_WORLD,"Fk     = %f\n", cell->Fk); CHKERRQ(ierr);
  ierr = PetscPrintf(MPI_COMM_WORLD,"Fk0    = %f\n", cell->Fk0); CHKERRQ(ierr);
  ierr = PetscPrintf(MPI_COMM_WORLD,"kclip  = %f\n", cell->kclip); CHKERRQ(ierr);
  ierr = PetscPrintf(MPI_COMM_WORLD,"ecm    = %f\n", cell->ecm); CHKERRQ(ierr);
  ierr = PetscPrintf(MPI_COMM_WORLD,"contactThres = %f\n", cell->contactThres); CHKERRQ(ierr);
  ierr = PetscPrintf(MPI_COMM_WORLD,"adhesionRadius = %f\n", cell->adhesionRadius); CHKERRQ(ierr);
  ierr = PetscPrintf(MPI_COMM_WORLD,"---------------\n", cell->ecm); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MyCellUpdateFluidFieldRHS"
PetscErrorCode MyCellUpdateFluidFieldRHS( DCell dcell, IIM iim, int ga, PetscReal t ) {
  PetscErrorCode ierr;
  PetscFunctionBegin;
  ierr = IIMSetForceComponents(iim, InterfacialForceAdhesion ); CHKERRQ(ierr);
  ierr = IIMSetForceContext(iim, dcell); CHKERRQ(ierr);
  MyCell mycell = (MyCell)dcell;
  mycell->sidx = iim->sidx;
  ierr = IIMUpdateRHS(iim, dcell->lsPlasmaMembrane, ga); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MyCellUpdateFluidFieldImplicitRHS"
PetscErrorCode MyCellUpdateFluidFieldImplicitRHS( DCell dcell, IIM iim, int ga, PetscReal t ) {
  PetscErrorCode ierr;
  PetscFunctionBegin;
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
  ierr = LevelSetWriteIrregularNodeList(dcell->lsPlasmaMembrane, ti); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MyCellWriteImplicit"
PetscErrorCode MyCellWriteImplicit( DCell dcell, int ti ) {
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
  PetscReal ***phi;
  iCoor p,q;
  Coor X;
  const int idh = 10; // samples per grid cell
  const PetscReal dh = 1. / idh;
  const Coor dH = ls->phi->d;
  const PetscReal dA = dh*dh*dh * dH.x*dH.y*dH.z;
  const PetscReal dx = ls->phi->d.x;
  PetscReal dist;
  int m, i, j, k;
  const int len = ArrayLength(ls->band);
  const iCoor *band = ArrayGetData(ls->band);
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = GridGetBounds(ls->phi,&p,&q); CHKERRQ(ierr);
  ierr = GridGet(ls->phi,&phi); CHKERRQ(ierr);
  cell->contactArea = 0;
  for (m = 0; m < len; ++m) {
    for (k = 0; k < idh; k++) {
      for (j = 0; j < idh; j++) {
        for (i = 0; i < idh; i++) {
          X.x = band[m].x + i * dh;
          X.y = band[m].y + j * dh;
          X.z = band[m].z + k * dh;
          ierr = GridInterpolate( lsGrooves->phi, X, &dist ); CHKERRQ(ierr);
          if( dist*dx > -cell->contactThres ) {
            cell->contactArea += LevelSetDiracDelta3D( phi, ls->phi->d, X ) * dA;
          } // if in contact
        } // i
      } // j
    } // k
  } // m in band
  ierr = PetscViewerASCIIPrintf(cell->contactareafile,"%e %e\n", t, cell->contactArea); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


void InterfacialForceGrooveAdhesion(IrregularNode *n, void *context )
{
  const MyCell c = (MyCell)context;
  const Coor dh = lsGrooves->phi->d;
  PetscReal dist;

  GridInterpolate( lsGrooves->phi, n->X, &dist );

  if( dist*dh.x > -c->contactThres ) {
    n->fa1 = c->Fa;
    n->fa2 = 0;
  } else {
    n->fa1 = 0;
    n->fa2 = 0;
  }

  PetscReal k = n->k;
  const PetscReal clip = c->kclip;
  k = k >  clip ?  clip : k;
  k = k < -clip ? -clip : k;

  n->F1 = c->scale * ( n->fa1 - c->Fk0 * k);
  n->F2 = c->scale * ( n->fa2 );
}
