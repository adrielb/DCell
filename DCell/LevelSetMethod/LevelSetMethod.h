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

//TODO: separate node coordinates from implementation dependent nodal quantities
typedef struct _IrregularNode IrregularNode;
struct _IrregularNode {
  iCoor pos;   //Node index in the grid
  int axis;       // 0: x-axis, 1: y-axis, 2: z-axis, -1: no axis
  VelFace shift;    //Irregnodes parallel to cell face U, V or W;
  int signCenter;  //sign of phi at pos
  int signFace;    //sign of phi at (pos-1)----|----pos
  int GAP; // TODO: hack, 4-byte gap so that bytes align in both 32 and 64-bit word boundaries
  // Coor o, n, s, r;
//  PetscReal ox, oy, oz; //Stencil intersection (similar to 'd' and 'axis') [may remove since redundant]
  PetscReal nx, ny, nz;//Normal direction
  PetscReal sx, sy, sz;//Tangential direction
  PetscReal rx, ry, rz;//Tangential direction (zero in 2D)
  PetscReal d;        //distance from central node
  PetscReal k;        //Curvature in 2D or 3D
  PetscReal k_nn, k_tt, k_nt;  // Principal Curvatures (k_nn used in 2D)
  PetscReal f1, f1_n, f1_t, f1_nn, f1_tt, f1_nt; // the normal force
  PetscReal f2, f2_n, f2_t, f2_nn, f2_tt, f2_nt; // the tangential force in local coor
  PetscReal f3, f3_n, f3_t, f3_nn, f3_tt, f3_nt;
  PetscReal ftx, fty, ftz;   // tangential force in global coor
  PetscReal fa1, fa2; // Force of adhesion (in normal/tangential direction)
  Coor X;   // X.x = x + ox + shift.x/2  (absolute position in space)
  Coor op;  //Ortho proj relative to node location
  int numNei;
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
inline PetscReal LevelSetDiracDelta2D( PetscReal **phi, const Coor dh, const Coor X );

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
