#ifndef FLUIDFIELD_PRIVATE_H_
#define FLUIDFIELD_PRIVATE_H_

//Log Events
PetscLogEvent EVENT_FluidFieldMaxVelocityMag;
PetscLogEvent EVENT_FluidField_DiscreteCompatibilityCondition;

//Private Functions
PetscErrorCode FluidFieldMatAssemble( FluidField f );
PetscErrorCode FluidField_MatAssemble( PetscReal mu, Array dbc, DA daV, Coor dH, Mat mat );
PetscErrorCode FluidFieldIntegrateStrainRate( DA daV, Vec vecV, DA daE, Vec vecE, PetscReal dh, PetscReal dt );
PetscErrorCode FluidFieldElasticDivergence( DA daE, Vec et, DA daV, Vec rhs, PetscReal dh );
PetscErrorCode FluidField_EnforceNoSlipBC( FluidField f );
PetscErrorCode FluidField_PressureBC( FluidField f );
PetscErrorCode FluidField_AppendDBC( Array dbc, MatStencil row );
PetscErrorCode FluidFieldRegisterEvents(  );

#endif /* FLUIDFIELD_PRIVATE_H_ */
