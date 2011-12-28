#include "DCell.h"

//TODO: delete this file?

#undef __FUNCT__
#define __FUNCT__ "RxnCreateNull"
PetscInt EVENT_RxnCreateNull;
PetscErrorCode RxnCreateNull( PetscInt num, Rxn *rxn )
{
  struct _Rxn *r;
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  PetscLogEventBegin(EVENT_RxnCreateNull,0,0,0,0);
//  PetscLogEventRegister(&EVENT_RxnCreateNull,"RxnCreateNull", 0);
  
  ierr = PetscNew(struct _Rxn, &r); CHKERRQ(ierr);  // Allocate space for reaction structure
  
  ierr = PetscMalloc(num*sizeof(PetscReal), &r->diffusion_coeff); CHKERRQ(ierr);
  ierr = PetscMalloc(num*sizeof(ReactionKineticsFunction), &r->R); CHKERRQ(ierr);
  ierr = PetscMalloc(num*sizeof(ReactionKineticsFunction), &r->dRdC); CHKERRQ(ierr);
  ierr = PetscMalloc(num*sizeof(PetscReal*), r->jac);; CHKERRQ(ierr);
  for( int i = 0; i < num; i++)
  {
    ierr = PetscMalloc(num*sizeof(PetscReal), ); CHKERRQ(ierr);
  }
  
  r->num = num;
  *rxn = r; // Return contructed object
  
  PetscLogEventEnd(EVENT_RxnCreateNull,0,0,0,0);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "RxnDestroy"
PetscInt EVENT_RxnDestroy;
PetscErrorCode RxnDestroy( Rxn rxn )
{
  PetscErrorCode ierr;
  
  
  PetscFunctionBegin;
  PetscLogEventBegin(EVENT_RxnDestroy,0,0,0,0);
//  PetscLogEventRegister(&EVENT_RxnDestroy,"RxnDestroy", 0);
  PetscFree(rxn->diffusion_coeff);
  PetscFree(rxn->R);
  PetscFree(rxn->jac);
  PetscFree(rxn);
  PetscLogEventEnd(EVENT_RxnDestroy,0,0,0,0);
  PetscFunctionReturn(0);
}

PetscReal Null_R( PetscReal *c)
{
  return 0;
}

PetscReal Null_dR_dC( PetscReal *c )
{
  return 0;
}