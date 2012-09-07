#ifndef LEVELSETMETHOD_H_
#define LEVELSETMETHOD_H_

#include "Grid.h"
#include "Common.h"

typedef struct _ParticleLS *ParticleLS;
typedef struct _LevelSet *LevelSet;

struct _LevelSet {
  Grid phi;  // the level set
  Grid phi0; //
  LevelSet psi;  // implicit step
  Grid tmp; // temporary buffer while resizing     // TODO: use a global set of work grids for all level sets
  Array band; // Narrow band coordinates [x,y,z]
  Array irregularNodes;
  PetscReal bandWidth; // +/- distance to the zero level set // TODO: bandWidth and PHI_INF are redundant constants! use one
  PetscReal PHI_INF;   // max distance outside narrow band
  PetscReal CFLthres;  // number of grid points ls moves before it triggers reinitialization
  PetscReal CFLcount;  // number of grid points ls has moved between reinitializations
  PetscInt AdvectThres;// number of times ls moves before it triggers reinitializations
  PetscInt AdvectCount;// number of times ls has moved between reinitializations
  PetscReal maxVel;
  ParticleLS pls;
  PetscErrorCode (*Advect)(LevelSet ls, Grid velgrid, PetscReal dt);
};

typedef struct _IrregularNode IrregularNode;
struct _IrregularNode {
  iCoor pos; // Node index in the grid
  Coor X;    // absolute position in space
  Coor op;   // Ortho proj relative to node location
  int sign;
};

/* Level Set */
PetscErrorCode LevelSetCreate(Coor dh, iCoor pos, iCoor size, LevelSet *levelset);
PetscErrorCode LevelSetDestroy(LevelSet ls);
PetscErrorCode LevelSetSetBandWidth(LevelSet ls, PetscReal bandwidth);
PetscErrorCode LevelSetDuplicate( LevelSet ls, LevelSet *copy);
PetscErrorCode LevelSetResize( LevelSet ls);
PetscErrorCode LevelSetUpdateIrregularNodeList( LevelSet ls );
PetscErrorCode LevelSetWriteIrregularNodeList( LevelSet ls, int idx );
PetscErrorCode LevelSetNormalDirection( LevelSet ls, Coor X, Coor *n );
inline PetscReal LevelSetDiracDelta2D( Grid phi, const Coor X );
inline PetscReal LevelSetDiracDelta3D( Grid phi, const Coor X );

PetscErrorCode LevelSetAdvect( LevelSet ls, int ga, PetscReal dt );
PetscErrorCode LevelSetAdvectSL(LevelSet ls, Grid velgrid, PetscReal dt);
PetscErrorCode LevelSetAdvectAndReinit(LevelSet ls, Grid velgrid, PetscReal dt);
PetscErrorCode LevelSetAdvectPLS(LevelSet ls, Grid velgrid, PetscReal dt);
PetscErrorCode LevelSetAdvectSLRK2HalfStep( LevelSet ls, Grid velgrid, PetscReal dt );
PetscErrorCode LevelSetAdvectSLRK2FullStep( LevelSet ls, Grid velgrid, PetscReal dt );
PetscErrorCode LevelSetAdvectImplicit( LevelSet ls, Grid velgrid, PetscReal dt );
PetscErrorCode LevelSetAdvectImplicitInit( LevelSet ls, PetscInt *n );
PetscErrorCode LevelSetAdvectImplicitRHS( LevelSet ls, int ga, PetscReal dt, PetscReal *g );
PetscErrorCode LevelSetAdvectImplicitUpdate( LevelSet ls, PetscReal lambda, PetscReal *dpsi );
PetscErrorCode LevelSetAdvectImplicitReinit( LevelSet ls, PetscReal dt );

PetscErrorCode LevelSetInitializeParticles( LevelSet ls );
PetscErrorCode ParticleLSDestroy( ParticleLS pls );
PetscErrorCode ParticleLSWriteParticles(ParticleLS pls, int t);
PetscErrorCode LevelSetAdvectPLSRK2HalfStep( LevelSet ls, Grid velgrid, PetscReal dt );
PetscErrorCode LevelSetAdvectPLSRK2FullStep( LevelSet ls, Grid velgrid, PetscReal dt );

PetscErrorCode LevelSetReinitialize( LevelSet ls );


/* Level Set Initializations */
PetscErrorCode LevelSetInitializeToCircle( Coor dh, Coor center, PetscReal radius, LevelSet *lset );
PetscErrorCode LevelSetInitializeToSphere( Coor dh, Coor center, PetscReal radius, LevelSet *lset );
PetscErrorCode LevelSetInitializeToStar2D( Coor dh, Coor center, PetscReal radius, PetscReal amp, PetscReal numPetals, LevelSet *lset );
PetscErrorCode LevelSetInitializeToStar3D( Coor dh, Coor center, PetscReal radius, PetscReal amp, PetscReal numPetals, LevelSet *lset );
PetscErrorCode LevelSetInitializeFromImage( LevelSet ls );
#endif /*LEVELSETMETHOD_H_*/
