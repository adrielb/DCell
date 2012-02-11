#ifndef DWORLD_H_
#define DWORLD_H_

#include "FluidField.h"

typedef struct _DWorld *DWorld;
typedef struct _DCell *DCell;

struct _DCell {
  int id;
  int type;
  LevelSet lsPlasmaMembrane;
  PetscErrorCode (*Destroy)(DCell this);
  PetscErrorCode (*Write)(DCell this, int t);
  PetscErrorCode (*UpdateFluidFieldRHS)(DCell this, IIM, int ga, PetscReal t );
  PetscErrorCode (*Advect)(DCell this, int ga, PetscReal dt );
  PetscErrorCode (*AdvectRK2HalfStep)(DCell this, int ga, PetscReal dt );
  PetscErrorCode (*AdvectRK2FullStep)(DCell this, int ga, PetscReal dt );
  PetscErrorCode (*AdvectImplicitInit)(DCell this, PetscInt *n );
  PetscErrorCode (*AdvectImplicitRHS)(DCell this, int ga, PetscReal dt, PetscReal *g );
  PetscErrorCode (*AdvectImplicitUpdate)(DCell this, PetscReal lambda, PetscReal *dpsi );
};

typedef struct _DCellsArray {
  Array dcells;
  Array sendIdx;
  Array sendRank;
  int totSend;

} *DCellsArray;

struct _DWorld {
  DCellsArray dcells; // cells this processor owns
  FluidField fluid;
  IIM iim;

  // Temporal stats
  PetscBool printStep;
  PetscReal maxVel;// current maximum velocity (magnitude)
  PetscReal CFL;   // current CFL constraint
  PetscReal dt;    // current time step
  PetscReal dtmax; // maximum time step
  PetscReal dtcfl; // maximum dt given CFL constraint
  PetscReal dtframe; // fixed time step for output animation data
  PetscInt  tiframe; // frame number
  PetscReal t;     // current time in units of time
  PetscReal tend;  // maximum simulation time
  PetscInt  ti;    // current time in iteration number
  PetscInt  timax; // maximum number of iterations
  PetscInt  writeInterval;

  PetscViewer temporalfile;

  PetscLogStage stageSimLoop;

  PetscErrorCode (*Simulate)( DWorld world );
};

// DWorld
PetscErrorCode DWorldCreate( FluidField fluid, DWorld *world );
PetscErrorCode DWorldDestroy( DWorld world );
PetscErrorCode DWorldAddDCell( DWorld world, void *dcell );
PetscErrorCode DWorldSimulate( DWorld world );
PetscErrorCode DWorldSetPrintStep( DWorld world, PetscBool printStep );
PetscErrorCode DWorldPrintStep( DWorld world );
PetscErrorCode DWorldWrite( DWorld world, int ti );
PetscErrorCode DWorldSetFrameInterval( DWorld w, PetscReal dtframe );
PetscErrorCode DWorldSetFromOptions( DWorld w );
PetscErrorCode DWorldSimulate_Euler( DWorld world );
PetscErrorCode DWorldSimulate_RK2( DWorld world );
PetscErrorCode DWorldSimulate_BFGS( DWorld world );

// DCells Array
PetscErrorCode DCellsArrayCreate(DCellsArray *dcellsarray);
PetscErrorCode DCellsArrayDestroy(DCellsArray dcells );
PetscErrorCode DCellsArrayAdd( DCellsArray dcells, DCell cell );
PetscErrorCode DCellsArrayAdvect( DCellsArray dcells, FluidField f, PetscReal dt );
PetscErrorCode DCellsArrayWrite(  DCellsArray dcells, int t);
PetscErrorCode DCellsArrayUpdateFluidFieldRHS( DCellsArray dcells, IIM iim, FluidField f, PetscReal t );
PetscErrorCode DCellsArrayAdvectRK2HalfStep( DCellsArray dcells, FluidField f, PetscReal dt );
PetscErrorCode DCellsArrayAdvectRK2FullStep( DCellsArray dcells, FluidField f, PetscReal dt );
PetscErrorCode DCellsArrayAdvectImplicit( DCellsArray dcells, FluidField f, PetscReal dt, PetscReal *p, PetscReal *g );
PetscErrorCode DCellsArrayAdvectImplicitInit( DCellsArray dcells, PetscInt *n );
PetscErrorCode DCellsArrayAdvectImplicitRHS( DCellsArray dcells, FluidField f, PetscReal dt, PetscReal *g );
PetscErrorCode DCellsArrayAdvectImplicitUpdate( DCellsArray dcells, PetscReal lambda, PetscReal *dpsi );

// DCell
PetscErrorCode DCellCreate( LevelSet lsPlasmaMembrane, DCell *dcell );
PetscErrorCode DCellSetup( LevelSet lsPlasmaMembrane, DCell cell );
PetscErrorCode DCellDestroy( DCell c );
PetscErrorCode DCellWrite( DCell dcell, int ti );
PetscErrorCode DCellAdvect( DCell dcell, int ga, PetscReal dt );
PetscErrorCode DCellAdvectRK2HalfStep( DCell dcell, int ga, PetscReal dt );
PetscErrorCode DCellAdvectRK2FullStep( DCell dcell, int ga, PetscReal dt );
PetscErrorCode DCellUpdateFluidFieldRHS( DCell this, IIM, int ga, PetscReal t );
PetscErrorCode DCellAdvectImplicitInit( DCell this, PetscInt *n );
PetscErrorCode DCellAdvectImplicitRHS( DCell this, int ga, PetscReal dt, PetscReal *g );
PetscErrorCode DCellAdvectImplicitUpdate( DCell this, PetscReal lambda, PetscReal *dpsi );
#endif /*DWORLD_H_*/
