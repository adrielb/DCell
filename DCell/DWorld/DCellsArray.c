#include "DWorld.h"

PetscErrorCode DCellsRegisterEvents( );

#undef __FUNCT__
#define __FUNCT__ "DCellsArrayCreate"
PetscErrorCode DCellsArrayCreate(DCellsArray *dcellsarray)
{
  DCellsArray cells;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscNew(struct _DCellsArray, &cells); CHKERRQ(ierr);
  ierr = ArrayCreate("dcell_array", sizeof(DCell), 10, &cells->dcells); CHKERRQ(ierr);
  *dcellsarray = cells;
  ierr = DCellsRegisterEvents(); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DCellsArrayDestroy"
PetscErrorCode DCellsArrayDestroy( DCellsArray dcells )
{
  int i;
  DCell dcell;
  PetscErrorCode ierr;
  for ( i = 0; i < ArrayLength(dcells->dcells); ++i) {
    ierr = ArrayGetP(dcells->dcells, i, &dcell);    CHKERRQ(ierr);
    ierr = dcell->Destroy(dcell); CHKERRQ(ierr);
  }
  ierr = ArrayDestroy(dcells->dcells);  CHKERRQ(ierr);
  ierr = PetscFree(dcells); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DCellsArrayAdd"
PetscErrorCode DCellsArrayAdd( DCellsArray dcells, DCell cell )
{
  DCell *newCell;
  PetscErrorCode ierr;
  ierr = ArrayAppend( dcells->dcells, &newCell); CHKERRQ(ierr);
  *newCell = cell;
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DCellsArrayAdvect"
static PetscLogEvent EVENT_DCellsArrayAdvect;
PetscErrorCode DCellsArrayAdvect(DCellsArray dcells, FluidField f, PetscReal dt )
{
  int i;
  DCell dcell;
  PetscErrorCode ierr;
  ierr = PetscLogEventBegin(EVENT_DCellsArrayAdvect,0,0,0,0); CHKERRQ(ierr);
  ierr = GAPutVec(f->vel,f->ga); CHKERRQ(ierr);
  for ( i = 0; i < ArrayLength(dcells->dcells); ++i) {
    ierr = ArrayGetP(dcells->dcells, i, &dcell); CHKERRQ(ierr);
    dcell->Advect(dcell, f->ga, dt);
  }
  ierr = PetscBarrier(0); CHKERRQ(ierr);
  ierr = PetscLogEventEnd(EVENT_DCellsArrayAdvect,0,0,0,0); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DCellsArrayAdvectRK2HalfStep"
static PetscLogEvent EVENT_DCellsArrayAdvectRK2HalfStep;
PetscErrorCode DCellsArrayAdvectRK2HalfStep(DCellsArray dcells, FluidField f, PetscReal dt )
{
  int i;
  DCell dcell;
  PetscErrorCode ierr;
  ierr = PetscLogEventBegin(EVENT_DCellsArrayAdvectRK2HalfStep,0,0,0,0); CHKERRQ(ierr);
  ierr = GAPutVec(f->vel,f->ga); CHKERRQ(ierr);
  for ( i = 0; i < ArrayLength(dcells->dcells); ++i) {
    ierr = ArrayGetP(dcells->dcells, i, &dcell); CHKERRQ(ierr);
    dcell->AdvectRK2HalfStep(dcell, f->ga, dt);
  }
  ierr = PetscBarrier(0); CHKERRQ(ierr);
  ierr = PetscLogEventEnd(EVENT_DCellsArrayAdvectRK2HalfStep,0,0,0,0); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DCellsArrayAdvectRK2FullStep"
static PetscLogEvent EVENT_DCellsArrayAdvectRK2FullStep;
PetscErrorCode DCellsArrayAdvectRK2FullStep(DCellsArray dcells, FluidField f, PetscReal dt )
{
  int i;
  DCell dcell;
  PetscErrorCode ierr;
  ierr = PetscLogEventBegin(EVENT_DCellsArrayAdvectRK2FullStep,0,0,0,0); CHKERRQ(ierr);
  ierr = GAPutVec(f->vel,f->ga); CHKERRQ(ierr);
  for ( i = 0; i < ArrayLength(dcells->dcells); ++i) {
    ierr = ArrayGetP(dcells->dcells, i, &dcell); CHKERRQ(ierr);
    dcell->AdvectRK2FullStep(dcell, f->ga, dt);
  }
  ierr = PetscBarrier(0); CHKERRQ(ierr);
  ierr = PetscLogEventEnd(EVENT_DCellsArrayAdvectRK2FullStep,0,0,0,0); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DCellsArrayWrite"
PetscErrorCode DCellsArrayWrite( DCellsArray dcells, int t)
{
  int i;
  DCell dcell;
  PetscErrorCode ierr;
  for ( i = 0; i < ArrayLength(dcells->dcells); ++i) {
    ierr = ArrayGetP(dcells->dcells, i, &dcell);    CHKERRQ(ierr);
    ierr = dcell->Write(dcell, t); CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DCellsArrayUpdateFluidFieldRHS"
static PetscLogEvent EVENT_DCellsArrayUpdateFluidFieldRHS;
PetscErrorCode DCellsArrayUpdateFluidFieldRHS( DCellsArray dcells, IIM iim, FluidField f, PetscReal t )
{
  int i;
  DCell dcell;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  // IIM update for each level set
  ierr = PetscLogEventBegin(EVENT_DCellsArrayUpdateFluidFieldRHS,0,0,0,0); CHKERRQ(ierr);
  ierr = GAPutVec(f->rhs,f->ga); CHKERRQ(ierr);
  for ( i = 0; i < ArrayLength(dcells->dcells); ++i) {
    ierr = ArrayGetP(dcells->dcells,i,&dcell); CHKERRQ(ierr);
    ierr = dcell->UpdateFluidFieldRHS(dcell, iim, f->ga, t ); CHKERRQ(ierr);
  }
  ierr = GAGetVec(f->ga,f->rhs); CHKERRQ(ierr);
  ierr = PetscLogEventEnd(EVENT_DCellsArrayUpdateFluidFieldRHS,0,0,0,0); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DCellsArrayAdvectImplicitInit"
PetscErrorCode DCellsArrayAdvectImplicitInit( DCellsArray dcells, int *n )
{
  int i;
  DCell dcell;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  *n = 0;
  for ( i = 0; i < ArrayLength(dcells->dcells); ++i) {
    ierr = ArrayGetP(dcells->dcells,i,&dcell); CHKERRQ(ierr);
    ierr = dcell->AdvectImplicitInit(dcell, n); CHKERRQ(ierr);
  }
  // TODO: sum up all sizes among procs
//  ierr = MPI_Allreduce(n,n,1,MPI_DOUBLE,MPI_SUM,comm); CHKERRQ(ierr);
  ierr = PetscBarrier(0); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DCellsArrayAdvectImplicitRHS"
PetscErrorCode DCellsArrayAdvectImplicitRHS( DCellsArray dcells, FluidField f, double dt, double *g0 )
{
  int i;
  DCell dcell;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = GAPutVec(f->vel,f->ga); CHKERRQ(ierr);
  for ( i = 0; i < ArrayLength(dcells->dcells); ++i) {
    ierr = ArrayGetP(dcells->dcells,i,&dcell); CHKERRQ(ierr);
    ierr = dcell->AdvectImplicitRHS(dcell, f->ga, dt, g0 ); CHKERRQ(ierr);
  }
  ierr = PetscBarrier(0); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DCellsArrayAdvectImplicitUpdate"
PetscErrorCode DCellsArrayAdvectImplicitUpdate( DCellsArray dcells, double lambda, double *d )
{
  int i;
  DCell dcell;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  for ( i = 0; i < ArrayLength(dcells->dcells); ++i) {
    ierr = ArrayGetP(dcells->dcells,i,&dcell); CHKERRQ(ierr);
    ierr = dcell->AdvectImplicitUpdate(dcell, lambda, d ); CHKERRQ(ierr);
  }
  ierr = PetscBarrier(0); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DCellsRegisterEvents"
static int EVENTS_registered = PETSC_FALSE;
PetscErrorCode DCellsRegisterEvents(  )
{
  PetscErrorCode ierr;
  if( EVENTS_registered )
    PetscFunctionReturn(0);

  ierr = PetscLogEventRegister("DCellsAdvect",0,&EVENT_DCellsArrayAdvect); CHKERRQ(ierr);
  ierr = PetscLogEventRegister("DCellsAdvect",0,&EVENT_DCellsArrayAdvectRK2HalfStep); CHKERRQ(ierr);
  ierr = PetscLogEventRegister("DCellsAdvect",0,&EVENT_DCellsArrayAdvectRK2FullStep); CHKERRQ(ierr);
  ierr = PetscLogEventRegister("DCellsUpdateRHS",0,&EVENT_DCellsArrayUpdateFluidFieldRHS); CHKERRQ(ierr);

  EVENTS_registered = PETSC_TRUE;
  PetscFunctionReturn(0);
}
