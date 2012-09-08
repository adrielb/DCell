//#ifndef IIM_PRIVATE_H_
//#define IIM_PRIVATE_H_
#include "ImmersedInterfaceMethod.h"

struct _IIM
{
  LocalCoor lc;// TODO: No need for this object, array of 12 elements is enough for stack
  LeastSq lsq; // TODO: No need for this object, array of 12 elements is enough for stack
  PetscBool is2D;
  PetscInt Np;
  PetscReal eps; // Radius around node considered a neighbor
  InterfacialForce F;  // Normal and tangential force functions
  void *context; // context passed to Interfacial force function
  Array idx,coor,val; //Interpolation indexes, coors and values
  Array debug;   // for iimtest.c
  PetscReal mu;  // Fluid viscosity
  Coor f;        // origin of fluid field coordinate system
  Coor df;       // grid widths of fluid field coordinate system
  Array iimIrregularNodes;
  SpatialIndex sidx; // Indexes irregular nodes in 3D bins
};

struct _LocalCoor {
  Array coor;
};

typedef struct IIMDebug {
  PetscReal h;
  PetscReal dd;
  PetscReal j, jp, jpp;
  PetscReal C;
  int dof;
  int sign;
  int axis;
  iCoor X;
} IIMDebug;

PetscLogEvent EVENT_IIMUpdateIrregularNodeGrid;
PetscLogEvent EVENT_IIMUpdateSurfaceQuantities;
PetscLogEvent EVENT_IIMUpdateSurfaceDerivatives;
PetscLogEvent EVENT_IIMUpdateRHS;
PetscLogEvent EVENT_IIMIrregularNodes;
PetscErrorCode IIMRegisterEvents( void );

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

PetscErrorCode IIMUpdateIrregularNodes( IIM iim, LevelSet ls );
PetscErrorCode IIMUpdateSurfaceQuantities( IIM iim, LevelSet ls );
PetscErrorCode IIMGridIntersections( IIM iim, const Grid g, const Coor p, const int axis );

void JumpPressure( IIMIrregularNode *n, Jump *j );
void JumpVelocity( PetscReal mu, IIMIrregularNode *n, Jump *j, int i );
PetscErrorCode IIMCorrection( IIM iim, iCoor x, Jump j, int axis, VelFace dof, int sign, PetscReal h, PetscReal dd, PetscReal mu );
PetscErrorCode IIMLaplaceCorrection( IIM iim, IIMIrregularNode *n, Jump j );
PetscErrorCode IIMPressureGradientCorrection( IIM iim, IIMIrregularNode *n, Jump j );
PetscErrorCode IIMVelocityGradientCorrection( IIM iim, IIMIrregularNode *n, Jump j );

PetscErrorCode IIMUpdateSurfaceQuantities_2D( IIM iim, LevelSet ls );
PetscErrorCode IIMUpdateSurfaceQuantities_3D( IIM iim, LevelSet ls );
PetscErrorCode IIMUpdateSurfaceDerivatives_2D( IIM iim, LevelSet ls );
PetscErrorCode IIMUpdateSurfaceDerivatives_3D( IIM iim, LevelSet ls );

void IIMRotation1st_2D( IIMIrregularNode *n, Jump *j );
void IIMRotation2nd_2D( IIMIrregularNode *n, Jump *j );
void IIMLocalToGlobal_1st( IIMIrregularNode *n, Jump *j);
void IIMLocalToGlobal_2nd( IIMIrregularNode *n, Jump *j);
void IIMGlobalToLocal_1st( IIMIrregularNode *n, Jump *j);
void IIMGlobalToLocal_2nd( IIMIrregularNode *n, Jump *j);

//#endif /* IIM_PRIVATE_H_ */
