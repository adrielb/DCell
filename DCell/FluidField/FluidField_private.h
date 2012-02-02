#ifndef FLUIDFIELD_PRIVATE_H_
#define FLUIDFIELD_PRIVATE_H_

//Log Events
PetscLogEvent EVENT_FluidFieldMaxVelocityMag;
PetscLogEvent EVENT_FluidField_DiscreteCompatibilityCondition;

//Private Functions
PetscErrorCode FluidFieldMatAssemble( FluidField f );
PetscErrorCode FluidField_MatAssemble( PetscReal mu, Array dbc, DM daV, Coor dH, Mat mat );
PetscErrorCode FluidFieldIntegrateStrainRate( DM daV, Vec vecV, DM daE, Vec vecE, PetscReal dh, PetscReal dt );
PetscErrorCode FluidFieldElasticDivergence( DM daE, Vec et, DM daV, Vec rhs, PetscReal dh );
PetscErrorCode FluidField_EnforceNoSlipBC( FluidField f );
PetscErrorCode FluidField_PressureBC( FluidField f );
PetscErrorCode FluidField_AppendDBC( Array dbc, MatStencil row );
PetscErrorCode FluidField_ComputeMatrix(DM dm,Vec x,Mat jac,Mat B,MatStructure *stflg);
PetscErrorCode FluidFieldRegisterEvents(  );

#endif /* FLUIDFIELD_PRIVATE_H_ */
