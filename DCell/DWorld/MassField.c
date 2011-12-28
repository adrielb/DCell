#include "MassField.h"

#undef __FUNCT__
#define __FUNCT__ "MassFieldJacobian"
PetscLogEvent EVENT_MassFieldJacobian;
PetscErrorCode MassFieldJacobian( MassField mf,  )
{
  PetscErrorCode ierr;
  
  
  PetscFunctionBegin;
  PetscLogEventBegin(EVENT_MassFieldJacobian,0,0,0,0);
//  PetscLogEventRegister("MassFieldJacobian", 0, &EVENT_MassFieldJacobian);
  
  
  PetscLogFlops( 0 );
  PetscLogEventEnd(EVENT_MassFieldJacobian,0,0,0,0);
  PetscFunctionReturn(0);
}