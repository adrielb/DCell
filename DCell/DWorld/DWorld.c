#include "DWorld.h"

PetscLogEvent EVENT_DWorldWrite;

#undef __FUNCT__
#define __FUNCT__ "DWorldCreate"
PetscErrorCode DWorldCreate( FluidField fluid, DWorld *world )
{
  struct _DWorld *w;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscNew(struct _DWorld, &w); CHKERRQ(ierr);
  
  char filename[PETSC_MAX_PATH_LEN];
  char wd[PETSC_MAX_PATH_LEN];
  ierr = DCellGetWorkingDirectory(wd); CHKERRQ(ierr);
  ierr = PetscSNPrintf(filename,PETSC_MAX_PATH_LEN,"%s/temporal.dat",wd); CHKERRQ(ierr);
  ierr = PetscViewerASCIIOpen(PETSC_COMM_SELF,filename,&w->temporalfile); CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(w->temporalfile,"iter t dt dtcfl tiframe CFL\n"); CHKERRQ(ierr);
  
  ierr = DCellsArrayCreate(&w->dcells); CHKERRQ(ierr);
  ierr = IIMCreate(!fluid->is3D,&fluid->mu,1.0,32,fluid->dh,&w->iim); CHKERRQ(ierr);

  w->fluid = fluid;
  w->timax = 1;
  w->tend = 1;
  w->CFL = 0.9;
  w->dt = 1;
  w->writeInterval = 1;
  w->dtframe = -1;
  w->tiframe = 0;
  w->printStep = PETSC_TRUE;
  w->Simulate = DWorldSimulate_BFGS;

  //BFGS
  int N = 100;
  MPI_Comm comm = PETSC_COMM_SELF;
  ierr = MatCreateSeqDense(comm,N,N,PETSC_NULL,&w->jac); CHKERRQ(ierr);
  ierr = VecCreateSeq(comm,N,&w->s); CHKERRQ(ierr);
  ierr = VecDuplicate(w->s,&w->y); CHKERRQ(ierr);
  ierr = VecDuplicate(w->s,&w->x0); CHKERRQ(ierr);
  ierr = VecDuplicate(w->s,&w->x1); CHKERRQ(ierr);

  ierr = PetscLogEventRegister("DWorldWrite", 0, &EVENT_DWorldWrite); CHKERRQ(ierr);
  ierr = PetscInfo(0, "Created DWorld\n"); CHKERRQ(ierr);

  *world = w;
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DWorldDestroy"
PetscErrorCode DWorldDestroy(DWorld world)
{
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  ierr = PetscViewerDestroy(&world->temporalfile); CHKERRQ(ierr);
  ierr = DCellsArrayDestroy(world->dcells); CHKERRQ(ierr);
  ierr = IIMDestroy(world->iim); CHKERRQ(ierr);
  ierr = FluidFieldDestroy(world->fluid); CHKERRQ(ierr);
  ierr = PetscFree(world); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DWorldAddDCell"
PetscErrorCode DWorldAddDCell( DWorld world, void *dcell )
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = DCellsArrayAdd( world->dcells, dcell ); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DWorldWrite"
PetscErrorCode DWorldWrite( DWorld world, int ti )
{
  char line[128];
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscLogEventBegin(EVENT_DWorldWrite,0,0,0,0); CHKERRQ(ierr);
  ierr = FluidFieldWrite( world->fluid,  ti); CHKERRQ(ierr);
  ierr = DCellsArrayWrite(world->dcells, ti); CHKERRQ(ierr);
  sprintf(line,"%d %e %e %e %d %e\n",
      world->ti,world->t,world->dt,world->dtcfl,world->tiframe, world->dt * world->maxVel / world->fluid->dh.x);
  ierr = PetscInfo1(0,"%s",line); CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(world->temporalfile,"%s",line); CHKERRQ(ierr);
  ierr = PetscLogEventEnd(EVENT_DWorldWrite,0,0,0,0); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DWorldSimulate"
PetscErrorCode DWorldSimulate(DWorld w)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = w->Simulate( w ); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode DWorldPrintStep( DWorld w )
{
  PetscErrorCode ierr;
  if( w->dtframe > 0 ) {
    ierr = PetscPrintf(PETSC_COMM_WORLD, "frame:\t %8d \n", w->tiframe-1);
  }
  ierr = PetscPrintf(PETSC_COMM_WORLD, "iter: \t %8d \t %8d \t %3.2f%%\n", w->ti, w->timax, (100.*w->ti)/w->timax);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "time: \t %.3e \t %.3e \t %3.2f%%\n", w->t, w->tend, (100.*w->t)/w->tend);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "dt:   \t %.3e \t %.3e \n", w->dt, w->dtcfl);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "CFL:  \t %f \t %f \n", w->dt*w->maxVel/w->fluid->dh.x, w->CFL);
  return 0;
}

PetscErrorCode DWorldSetPrintStep( DWorld w, PetscBool printStep )
{
  w->printStep = printStep;
  return 0;
}

PetscErrorCode DWorldSetWriteInterval( DWorld w, PetscInt interval )
{
  w->writeInterval = interval;
  return 0;
}

PetscErrorCode DWorldSetFrameInterval( DWorld w, PetscReal dtframe )
{
  w->dtframe = dtframe;
  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "DWorldSetFromOptions"
PetscErrorCode DWorldSetFromOptions( DWorld w )
{
  MPI_Comm comm = w->fluid->comm;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscOptionsGetReal(0,"-CFL",&w->CFL,0); CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(0,"-dtmax",&w->dtmax,0); CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(0,"-dtframe",&w->dtframe,0); CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(0,"-tend",&w->tend,0); CHKERRQ(ierr);
  ierr = PetscOptionsGetInt (0,"-timax",&w->timax,0); CHKERRQ(ierr);
  ierr = PetscOptionsGetInt (0,"-writeInterval",&w->writeInterval,0); CHKERRQ(ierr);

  ierr = PetscPrintf(comm, "DWorld Options\n"); CHKERRQ(ierr);
  ierr = PetscPrintf(comm, "CFL     = %e\n", w->CFL); CHKERRQ(ierr);
  ierr = PetscPrintf(comm, "dtmax   = %e\n", w->dtmax); CHKERRQ(ierr);
  ierr = PetscPrintf(comm, "dtframe = %e\n", w->dtframe); CHKERRQ(ierr);
  ierr = PetscPrintf(comm, "tend    = %e\n", w->tend); CHKERRQ(ierr);
  ierr = PetscPrintf(comm, "timax   = %d\n", w->timax); CHKERRQ(ierr);
  ierr = PetscPrintf(comm, "writeInterval   = %d\n", w->writeInterval); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
