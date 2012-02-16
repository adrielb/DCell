#ifndef FIELD_H_
#define FIELD_H_

#include "ImmersedInterfaceMethod.h"

typedef struct {
  PetscReal u,v,w,p;
} VelocityVector;

typedef struct {
  PetscReal xx,xy,xz,
               yy,yz,
                  zz;
} StrainTensor;

typedef struct _FluidField {
  // DA for velocity [u v w p]
  // [mat] {vel} = {rhs}
  DM daV;
  Mat mat;
  KSP ksp;
  Vec rhs;
  Vec vel;  // velocity at t = n+1
  Vec vel0; // velocity at t = n
  Array dirichletBC;
  int ga;

  /*
   * DA for strain tensor
   *  [exx exy exz]
   *  [exy eyy eyz]
   *  [exz eyz ezz]
   */
  DM daE;
  Vec E;

  // DA for single component buffer
  DM daB;
  Vec buf;

  PetscReal mu;
  Coor dh;   // Discretization size
  iCoor dims;
  PetscBool is3D;
  MPI_Comm comm;
  Grid mask;
} *FluidField;

//Public Functions
PetscErrorCode FluidFieldCreate(MPI_Comm comm, FluidField *fluid);
PetscErrorCode FluidFieldDestroy( FluidField f );
PetscErrorCode FluidFieldSetDims( FluidField f, iCoor dims );
PetscErrorCode FluidFieldSetDx(   FluidField f, PetscReal dx );
PetscErrorCode FluidFieldSetMask( FluidField f, Grid mask );
PetscErrorCode FluidFieldSetViscosity( FluidField f, PetscReal mu );
PetscErrorCode FluidFieldSetup( FluidField f );
PetscErrorCode FluidFieldWrite( FluidField f, int t);
PetscErrorCode FluidFieldSolve( FluidField f );
PetscErrorCode FluidFieldMaskDomain(FluidField f);
PetscErrorCode FluidFieldMaxVelocityMag( FluidField f, PetscReal *maxVel );

#endif /* FIELD_H_ */
