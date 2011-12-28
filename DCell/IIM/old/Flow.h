#ifndef FLOW_H_
#define FLOW_H_

#include "../Common/Main.h"
#include "petscdmmg.h"
#include "petscmg.h"

PetscInt EVENT_ComputeRHS;
PetscErrorCode ComputeRHS( DMMG dmmg, Vec b);

PetscInt EVENT_ComputeJacobian;
PetscErrorCode ComputeJacobian( DMMG dmmg,Mat jac,Mat B );




#endif /*FLOW_H_*/
