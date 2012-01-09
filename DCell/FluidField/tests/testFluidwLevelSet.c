#include "FluidField.h"

void InterfacialForceAdhesion(IrregularNode *n, void *context );
PetscErrorCode DCellsAdvect(FluidField f, IIM iim, Array lsets, PetscReal dt,
    int t);
PetscErrorCode DCellsWrite(Array lsets, int t);
PetscErrorCode DCellsDestroy(Array lsets);

typedef struct _Context {
  PetscReal K;  // magnitude of surface tension
  PetscReal fa; // magnitude of adhesion
  PetscReal scale;
  PetscReal contactThres;
  PetscReal contactArea;
  Coor dh;
  PetscViewer file;
} Context;

PetscErrorCode UpdateContext( Context *c, LevelSet ls, PetscReal t, PetscReal dt, PetscReal dtcfl )
{
  PetscErrorCode ierr;
  int len = ArrayLength(ls->irregularNodes);
  IrregularNode *nodes = ArrayGetData(ls->irregularNodes);

  c->contactArea = 0;
  for (int i = 0; i < len; ++i) {
    if( nodes[i].shift.x+nodes[i].shift.y != 0 ) continue; // if face node, skip it; (must be a cell-center node)
    if( nodes[i].X.y < c->contactThres )   c->contactArea += c->dh.x;
  }

  ierr = PetscViewerASCIIPrintf(c->file,"%e %e %e %e\n",t,c->contactArea,dt,dtcfl); CHKERRQ(ierr);

  return 0;
}

void InterfacialForceAdhesion(IrregularNode *n, void *context )
{
//  const Coor g = {0,-1,0}; //direction of adhesion normal to surface
  const Context *c = (Context*)context;

  if (n->X.y < c->contactThres) {
    n->fa1 = c->fa/c->dh.x; //direction of adhesion normal to membrane
    n->fa2 = 0;
//    n->fa1 = ( n->nx * g.x + n->ny * g.y)*c->fa/c->dh.x;
//    n->fa2 = (-n->ny * g.x + n->nx * g.y)*c->fa/c->dh.x;
  }

  n->f1 = c->scale*(n->fa1 - c->K * n->k);
  n->f2 = c->scale*(n->fa2);
}

int main(int argc, char **args) {
  PetscErrorCode ierr;
  ierr = DCellInit(); CHKERRQ(ierr);

  PetscReal d1 = 32, dx = 1./d1;
  iCoor size = {d1,d1,0};

  FluidField fluid;
  ierr = FluidFieldCreate(PETSC_COMM_WORLD, &fluid);  CHKERRQ(ierr);
  ierr = FluidFieldSetDims(fluid,size); CHKERRQ(ierr);
  ierr = FluidFieldSetDx(fluid,dx); CHKERRQ(ierr);
  ierr = FluidFieldSetup(fluid); CHKERRQ(ierr);

  IIM iim;
  ierr = IIMCreate(!fluid->is3D,&fluid->mu,1.5,32,fluid->dh,&iim); CHKERRQ(ierr);
  ierr = IIMSetForceComponents(iim,InterfacialForceAdhesion ); CHKERRQ(ierr);


  Array lsets;
  ierr = ArrayCreate("levelsets", sizeof(LevelSet), 1, &lsets);  CHKERRQ(ierr);

  int rank;
  MPI_Comm_rank(fluid->comm, &rank);
  if( rank == 0 ) {
    LevelSet *ls;
    ierr = ArrayAppend(lsets,&ls);  CHKERRQ(ierr);
//    ierr = LevelSetInitializeToStar2D(f->dh, (Coor) {.5,.5,0}, .4, .1, 8, ls); CHKERRQ(ierr);
//    ierr = ArrayWrite((*ls)->band,"band",0); CHKERRQ(ierr);
//    ierr = LevelSetInitializeToCircle(f->dh, (Coor) {.5,.4,0}, 0.35, ls); CHKERRQ(ierr);
    ierr = LevelSetInitializeToCircle(fluid->dh, (Coor) {0.5,0.5,0}, 0.4, ls); CHKERRQ(ierr);
  }

  PetscReal maxVel;
  PetscReal CFL = .1;
  PetscReal dt=0;
  PetscReal dtmax = 1e-3;
  PetscReal dtcfl = 0;
  PetscReal tend = 2e0;
  PetscInt  tmax = 2000;

  Context context;
  context.dh = fluid->dh;
  context.K = 10;
  context.fa = 100;
  context.scale = 1e-6;
  context.contactThres = 1.9;
  ierr = IIMSetForceContext(iim, &context); CHKERRQ(ierr);

  char filename[PETSC_MAX_PATH_LEN];
  char wd[PETSC_MAX_PATH_LEN];
  ierr = DCellGetWorkingDirectory(wd); CHKERRQ(ierr);
  ierr = PetscSNPrintf(filename,PETSC_MAX_PATH_LEN,"%s/temporal.dat",wd); CHKERRQ(ierr);
  ierr = PetscViewerASCIIOpen(PETSC_COMM_SELF,filename,&context.file); CHKERRQ(ierr);

  int t;
  for ( t = 0; t*dt < tend && t < tmax; ++t) {

    for (int i = 0; i < ArrayLength(lsets); ++i) {
      LevelSet ls;
      ierr = ArrayGetP(lsets, i, &ls);    CHKERRQ(ierr);
      ierr = UpdateContext(&context,ls, t*dt, dt, dtcfl); CHKERRQ(ierr);
    }

    ierr = FluidFieldSolve(fluid, iim, lsets, t); CHKERRQ(ierr);
    ierr = FluidFieldWrite(fluid, t); CHKERRQ(ierr);
    ierr = DCellsWrite(lsets,t); CHKERRQ(ierr);

    ierr = FluidFieldMaxVelocityMag(fluid,&maxVel); CHKERRQ(ierr);
    dtcfl = CFL * fluid->dh.x / maxVel;
    dt = PetscMin(dtmax,dtcfl);
    ierr = PetscPrintf(PETSC_COMM_WORLD, "time: \t %e \t %e \t %1.2f%%\n", t*dt, tend, 100*t*dt/tend);
    ierr = PetscPrintf(PETSC_COMM_WORLD, "iter: \t %d \t\t\t %1.1f \n", t, tend/dt);
    ierr = PetscPrintf(PETSC_COMM_WORLD, "dt:   \t %e \t %e \n", dt, dtcfl);
    ierr = PetscPrintf(PETSC_COMM_WORLD, "CFL:  \t %f \t %f \n", dt*maxVel/fluid->dh.x, CFL);

    ierr = DCellsAdvect(fluid, iim, lsets, dt, t+1); CHKERRQ(ierr);
  }
//  ierr = DCellsWrite(lsets,t+1); CHKERRQ(ierr);

  ierr = PetscViewerDestroy(context.file); CHKERRQ(ierr);

  ierr = IIMDestroy(iim); CHKERRQ(ierr);
  ierr = DCellsDestroy(lsets); CHKERRQ(ierr);
  ierr = FluidFieldDestroy(fluid); CHKERRQ(ierr);
  ierr = DCellFinalize(); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode DCellsAdvect(FluidField f, IIM iim, Array lsets, PetscReal dt,int t)
{
  LevelSet ls;
  PetscErrorCode ierr;

  ierr = GAPutVec(f->vel,f->ga); CHKERRQ(ierr);
  for (int i = 0; i < ArrayLength(lsets); ++i) {
    ierr = ArrayGetP(lsets, i, &ls);    CHKERRQ(ierr);
    ierr = LevelSetAdvect(ls, f->ga, dt);    CHKERRQ(ierr);
  }

  PetscBarrier(0);
  PetscFunctionReturn(0);
}

PetscErrorCode DCellsWrite(Array lsets, int t)
{
  LevelSet ls;
  PetscErrorCode ierr;
  for (int i = 0; i < ArrayLength(lsets); ++i) {
    ierr = ArrayGetP(lsets, i, &ls);    CHKERRQ(ierr);
//    ierr = GridWrite(ls->phi, t);    CHKERRQ(ierr);
    ierr = LevelSetWriteIrregularNodeList(ls->irregularNodes, t);    CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

PetscErrorCode DCellsDestroy(Array lsets)
{
  LevelSet ls;
  PetscErrorCode ierr;
  for (int i = 0; i < ArrayLength(lsets); ++i) {
    ierr = ArrayGetP(lsets, i, &ls);    CHKERRQ(ierr);
    ierr = LevelSetDestroy(ls);    CHKERRQ(ierr);
  }
  ierr = ArrayDestroy(lsets);  CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
