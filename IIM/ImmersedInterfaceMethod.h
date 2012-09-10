#ifndef IMMERSEDINTERFACEMETHOD_H_
#define IMMERSEDINTERFACEMETHOD_H_

#include "LevelSetMethod.h"

typedef struct _IIMIrregularNode IIMIrregularNode;
struct _IIMIrregularNode {
  SpatialItem spatialitem;
  PetscReal nx, ny, nz; // Normal direction
  PetscReal sx, sy, sz; // Tangential direction
  PetscReal rx, ry, rz; // Tangential direction (zero in 2D)
  PetscReal k;          // Curvature in 2D or 3D
  PetscReal k_nn, k_tt, k_nt;  // Principal Curvatures (k_nn used in 2D)
  PetscReal F1, f1, f1_n, f1_t, f1_nn, f1_tt, f1_nt; // the normal force
  PetscReal F2, f2, f2_n, f2_t, f2_nn, f2_tt, f2_nt; // the tangential force in local coor
  PetscReal F3, f3, f3_n, f3_t, f3_nn, f3_tt, f3_nt;
  PetscReal ftx, fty, ftz;   // tangential force in global coor
  PetscReal fa1, fa2; // Force of adhesion (in normal/tangential direction)
  Coor X;   // absolute position in space
  int numNei;
};

typedef void (*InterfacialForce)( IIMIrregularNode *n, void *context);
void InterfacialForceSurfaceTension( IIMIrregularNode *n, void *context );
void InterfacialForceElastic( IIMIrregularNode *n, void *context );
void InterfacialForceSimple( IIMIrregularNode *n, void *context );


typedef struct _IIM *IIM;
PetscErrorCode IIMCreate( PetscBool is2D, IIM *iim );
PetscErrorCode IIMDestroy( IIM iim );
PetscErrorCode IIMUpdateRHS( IIM iim, LevelSet ls, int ga );
PetscErrorCode IIMCorrectVelocity( IIM iim, const Coor X, Coor *vel );
PetscErrorCode IIMWriteIrregularNodeList( IIM iim, char *name, int idx );

PetscErrorCode IIMSetForceComponents( IIM iim, InterfacialForce F );
PetscErrorCode IIMSetForceContext( IIM iim, void *context );
PetscErrorCode IIMSetViscosity( IIM iim, PetscReal mu );
PetscErrorCode IIMSetEps( IIM iim, PetscReal eps );
PetscErrorCode IIMSetNp( IIM iim, int Np );
PetscErrorCode IIMSetFluidCoordinateSystem( IIM iim, Coor origin, Coor df );


typedef struct _LocalCoor *LocalCoor;
PetscErrorCode LocalCoorCreate( LocalCoor *lc );
PetscErrorCode LocalCoorDestroy( LocalCoor lc );
void LocalCoorSetLength( LocalCoor lc, PetscInt len);
void LocalCoorSolve( LocalCoor lc, Coor dh, IIMIrregularNode *N, IIMIrregularNode *nodes[] );
void LocalCoorGetVecs( LocalCoor lc, PetscReal **s, PetscReal **n, PetscReal **r );
void LocalCoor3DTangential( IIMIrregularNode *n );

#endif /*IMMERSEDINTERFACEMETHOD_H_*/
