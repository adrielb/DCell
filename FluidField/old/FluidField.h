#ifndef FLUIDFIELD_H_
#define FLUIDFIELD_H_

#include "ImmersedInterfaceMethod.h"

typedef struct _FluidField *FluidField;
typedef struct _FluidField {
  DA da;
  Vec source;    // Source-sink distribution 
  Mat matP;      // Poisson discretization 
  KSP kspP  ;    // Poisson Solver
  MatNullSpace nullspace; // The constant null space due to nueman bc for pressure
  Vec px,py,pz,p;// Pressure and pressure gradients
  Vec div;       // Divergence of velocity
  Vec u, v, w;   // Velocity vectors
  int gaU, gaV, gaW; // GlobalArray pointers for velocity
  Mat matU, matV, matW; // Velocity poisson discretization
  KSP kspU;      // Velocity solver
  PetscTruth is3D; // 2D or 3D simulation
  Coor l;        // Domain dimensions in units of length
  Coor d;        // Grid spacing
  PetscReal mu;  // Viscosity
  PetscReal rho; // Density
};

PetscErrorCode FluidFieldCreate(FluidField *fluid);
PetscErrorCode FluidFieldDestroy(FluidField f);
PetscErrorCode FluidFieldAssembleStaggeredGrid_2D( PetscReal dx, FluidField f );
PetscErrorCode FluidFieldAssembleStaggeredGrid_3D( PetscReal dx, FluidField f );
PetscErrorCode FluidFieldStep(FluidField f, IIM iim, int lsLEN, LevelSet *ls);
PetscErrorCode FluidFieldGetSource( FluidField f, Vec *source );

PetscErrorCode FluidFieldGradient( FluidField f );
PetscErrorCode FluidFieldGradient_2D(DALocalInfo info, Coor d, PetscReal mu, PetscReal rho, Vec P, Vec S, Vec PX, Vec PY );
PetscErrorCode FluidFieldGradient_3D(DALocalInfo info, Coor d, PetscReal mu, PetscReal rho, Vec P, Vec S, Vec PX, Vec PY, Vec PZ );

PetscErrorCode FluidFieldDivergence( FluidField f);
PetscErrorCode FluidFieldDivergence_2D( DALocalInfo info, Coor d, Vec US, Vec VS, Vec DIV );
PetscErrorCode FluidFieldDivergence_3D( DALocalInfo info, Coor d, Vec US, Vec VS, Vec WS, Vec DIV );

void AssembleLaplacian_BoundaryCheck2D( Mat mat, MatStencil row, int* lens, 
    PetscReal dx2, PetscReal dy2 );


#endif /*FLUIDFIELD_H_*/
