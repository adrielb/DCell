#ifndef FIBERFIELD_PRIVATE_H_
#define FIBERFIELD_PRIVATE_H_

PetscErrorCode FiberField_Init( FiberField field );
PetscErrorCode FiberField_ToVec( FiberField field );
PetscErrorCode FiberField_Step( FiberField field );
PetscErrorCode FiberField_Solve( FiberField f, PetscReal ti, Vec x );
PetscErrorCode FiberField_BodyCollisions( FiberField f );
PetscErrorCode FiberField_AssembleDisplacementMatrix( FiberField field );
PetscErrorCode FiberField_AssembleForceMatrix( FiberField field );
PetscErrorCode FiberField_FromVec( FiberField field );

#endif /* FIBERFIELD_PRIVATE_H_ */
