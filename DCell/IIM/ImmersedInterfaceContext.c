#include "ImmersedInterfaceMethod.h"
#include "IIM_private.h"

#undef __FUNCT__
#define __FUNCT__ "IIMCreate"
PetscErrorCode IIMCreate( PetscBool is2D, PetscReal *mu, PetscReal eps, int Np, Coor dh, IIM *iim )
{
  PetscErrorCode ierr;
  IIM i;
  
  PetscFunctionBegin;
  
  ierr = PetscNew(struct _IIM, &i); CHKERRQ(ierr);
  ierr = LocalCoorCreate(Np, &i->lc); CHKERRQ(ierr);
  ierr = LeastSqCreate(Np, is2D, &i->lsq); CHKERRQ(ierr);
  
  int est = 1e5;
  ierr = ArrayCreate( "iim_idx", sizeof(int*), est, &i->idx); CHKERRQ(ierr);
  ierr = ArrayCreate( "iim_coor", 4*sizeof(int), est, &i->coor); CHKERRQ(ierr); // [ x y z d ]
  ierr = ArrayCreate( "iim_val", sizeof(PetscReal), est, &i->val); CHKERRQ(ierr);
  ierr = ArrayCreate( "iim_grid",sizeof(GridPoint), 1, &i->irregularNodeGrid); CHKERRQ(ierr);
  //ierr = PetscOptionsGetReal(0,"-iim_eps",&t->eps,0); CHKERRQ(ierr);
//  ierr = PetscOptionsGetInt(0,"-iim_Np",&t->Np,0); CHKERRQ(ierr);
  i->dh = dh;
  i->mu = mu;
  i->eps = eps;
  i->Np = Np; 
  i->F = InterfacialForceSurfaceTension;
  
  *iim = i;

  ierr = IIMRegisterEvents(); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "IIMDestroy"
PetscErrorCode IIMDestroy( IIM iim )
{
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  ierr = ArrayDestroy(iim->coor); CHKERRQ(ierr);
  ierr = ArrayDestroy(iim->val); CHKERRQ(ierr);
  ierr = ArrayDestroy(iim->idx); CHKERRQ(ierr);
  ierr = ArrayDestroy(iim->irregularNodeGrid); CHKERRQ(ierr);
  ierr = LocalCoorDestroy( iim->lc ); CHKERRQ(ierr);
  ierr = LeastSqDestroy( iim->lsq ); CHKERRQ(ierr);
  ierr = PetscFree(iim); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode IIMSetForceComponents(IIM iim, InterfacialForce F )
{
  iim->F = F;
  return 0;
}

PetscErrorCode IIMSetForceContext(IIM iim, void *context )
{
  iim->context = context;
  return 0;
}

PetscErrorCode IIMSetNp( IIM iim, int Np )
{
  SETERRQ(PETSC_COMM_SELF,0,"NOT IMP");
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "IIMRegisterEvents"
static PetscBool EVENTS_registered = PETSC_FALSE;
PetscErrorCode IIMRegisterEvents()
{
  PetscErrorCode ierr;
  PetscFunctionBegin;
  if( EVENTS_registered )
    PetscFunctionReturn(0);

  ierr = PetscLogEventRegister("IIMUpdateRHS",0,&EVENT_IIMUpdateRHS); CHKERRQ(ierr);
  ierr = PetscLogEventRegister("IIMIrregGrid",0,&EVENT_IIMUpdateIrregularNodeGrid); CHKERRQ(ierr);
  ierr = PetscLogEventRegister("IIMSurfQuant",0,&EVENT_IIMUpdateSurfaceQuantities); CHKERRQ(ierr);
  ierr = PetscLogEventRegister("IIMSurfDeriv",0,&EVENT_IIMUpdateSurfaceDerivatives); CHKERRQ(ierr);

  EVENTS_registered = PETSC_TRUE;
  PetscFunctionReturn(0);
}
