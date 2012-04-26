#ifndef IMMERSEDINTERFACEMETHOD_H_
#define IMMERSEDINTERFACEMETHOD_H_

#include "LevelSetMethod.h"

typedef struct _LocalCoor *LocalCoor;
PetscErrorCode LocalCoorCreate( int Np, LocalCoor *lc );
PetscErrorCode LocalCoorDestroy( LocalCoor lc );
void LocalCoorSetLength( LocalCoor lc, PetscInt len);
void LocalCoor2DNormal( LocalCoor lc, LevelSet ls, IrregularNode *n );
void LocalCoor2DSolveStencil( LocalCoor lc, IrregularNode *n );
void LocalCoor2DGetVecs( LocalCoor lc, PetscReal **eta, PetscReal **xi );
void LocalCoor2DSolve( LocalCoor lc, Coor dh, IrregularNode *n);
void LocalCoor3DSolve( LocalCoor lc, Coor dh, IrregularNode *n );
void LocalCoor2DToArcLength( LocalCoor lc, IrregularNode *n,
    PetscReal n1, PetscReal n2, int i, PetscReal *s );
void LocalCoor3DGetVecs( LocalCoor lc, PetscReal **s, PetscReal **r, PetscReal **n );
void LocalCoor3DTangential( IrregularNode *n );


typedef void (*InterfacialForce)(IrregularNode *n, void *context);
void InterfacialForceSurfaceTension( IrregularNode *n, void *context );
void InterfacialForceElastic( IrregularNode *n, void *context );
void InterfacialForceSimple( IrregularNode *n, void *context );

typedef struct _IIM
{
  LocalCoor lc;// TODO: No need for this object, array of 12 elements is enough for stack
  LeastSq lsq; // TODO: No need for this object, array of 12 elements is enough for stack
  PetscInt Np;
  PetscReal eps; // Radius around node considered a neighbor
  PetscReal mu; // Fluid viscosity
//  iCoor dim;     // Dimension of FluidField
  InterfacialForce F;  //Normal and tangential force functions
  void *context; // context passed to Interfacial force function
  Array idx,coor,val; //Interpolation indexes, coors and values
  Coor dh;       // Grid widths
  SpatialIndex sidx; // Indexes irregular nodes in 3D bins

  iCoor p; // position of grid in space;
  iCoor n; // size of grid [H,W,L]
  Array irregularNodeGrid; // converts (x,y,z) to index in irregNodes array (TODO: implement this as a quad/octree or HRLE)
} *IIM;

typedef struct Jump
{
  PetscReal j; // Jump in the solution 
  
  //Local coordinate system jumps
  PetscReal e, n, t; 
  PetscReal ne, nn, nt, tt, te, ee;
  
  //Finite difference coordinate system jumps
  PetscReal x, y, z; // jump in first order derivatives
  PetscReal xx, yy, zz, xy, xz, yz; // jump in second order derivatives
} Jump;

typedef void (*JumpCondition)( PetscReal mu, IrregularNode *n, Jump *j, PetscInt axis_index);
void JumpPressure( IrregularNode *n, Jump *j );
void JumpVelocity(PetscReal mu, IrregularNode *n, Jump *j, int i );

PetscErrorCode IIMCreate( PetscBool is2D, int Np, Coor dh, IIM *iim );
PetscErrorCode IIMDestroy( IIM iim );
PetscErrorCode IIMSetForceComponents(IIM iim, InterfacialForce F );
PetscErrorCode IIMSetForceContext(IIM iim, void *context);
PetscErrorCode IIMSetViscosity( IIM iim, PetscReal mu );
PetscErrorCode IIMSetEps( IIM iim, PetscReal eps );
PetscErrorCode IIMUpdateRHS( IIM iim, LevelSet ls, int ga );

//Converts (x,y,z) to idx for spatial indexing
typedef struct {
  int idx;  // index of first node in irregularNode array at (x,y,z)
  int len;  // number of irregular nodes at (x,y,z)
} GridPoint;

PetscErrorCode IIMCorrection( IIM iim, iCoor x, Jump j, int axis, VelFace dof, int sign, PetscReal h, PetscReal dd, PetscReal mu );
PetscErrorCode IIMLaplaceCorrection( IIM iim, IrregularNode *n, Jump j );
PetscErrorCode IIMPressureGradientCorrection( IIM iim, IrregularNode *n, Jump j );
PetscErrorCode IIMVelocityGradientCorrection( IIM iim, IrregularNode *n, Jump j );
PetscErrorCode IIMUpdateSurfaceQuantities( IIM iim, LevelSet ls );
PetscErrorCode IIMUpdateSurfaceIndex_2D( IIM iim, LevelSet ls );
PetscErrorCode IIMUpdateSurfaceIndex_3D( IIM iim, LevelSet ls );
PetscErrorCode IIMUpdateSurfaceQuantities_2D( IIM iim, LevelSet ls );
PetscErrorCode IIMUpdateSurfaceQuantities_3D( IIM iim, LevelSet ls );
PetscErrorCode IIMUpdateSurfaceDerivatives_2D( IIM iim, LevelSet ls );
PetscErrorCode IIMUpdateSurfaceDerivatives_3D( IIM iim, LevelSet ls );

void IIMLocalToGlobal_1st( IrregularNode *n, Jump *j);
void IIMLocalToGlobal_2nd( IrregularNode *n, Jump *j);
void IIMGlobalToLocal_1st( IrregularNode *n, Jump *j);
void IIMGlobalToLocal_2nd( IrregularNode *n, Jump *j);

#endif /*IMMERSEDINTERFACEMETHOD_H_*/
