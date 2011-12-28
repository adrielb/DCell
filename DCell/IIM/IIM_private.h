#ifndef IIM_PRIVATE_H_
#define IIM_PRIVATE_H_

PetscLogEvent EVENT_IIMUpdateIrregularNodeGrid;
PetscLogEvent EVENT_IIMUpdateSurfaceQuantities;
PetscLogEvent EVENT_IIMUpdateSurfaceDerivatives;
PetscLogEvent EVENT_IIMUpdateRHS;

PetscErrorCode IIMRegisterEvents();

#endif /* IIM_PRIVATE_H_ */
