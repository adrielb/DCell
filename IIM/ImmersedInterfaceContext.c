#include "ImmersedInterfaceMethod.h"
#include "IIM_private.h"

#undef __FUNCT__
#define __FUNCT__ "IIMCreate"
PetscErrorCode IIMCreate( PetscBool is2D, int Np, Coor dh, IIM *iim )
{
  PetscErrorCode ierr;
  IIM i;
  
  PetscFunctionBegin;
  
  ierr = PetscNew(struct _IIM, &i); CHKERRQ(ierr);

  ierr = ArrayCreate( "iim_idx", sizeof(int*), &i->idx); CHKERRQ(ierr);
  ierr = ArrayCreate( "iim_coor", 4*sizeof(int), &i->coor); CHKERRQ(ierr); // [ x y z d ]
  ierr = ArrayCreate( "iim_val", sizeof(PetscReal), &i->val); CHKERRQ(ierr);

  i->dh = dh;
  i->mu = 1.0;
  i->eps = 1.1;
  i->Np = Np;
  i->F = InterfacialForceSurfaceTension;
  
  ierr = PetscOptionsGetReal(0, "-iim_eps", &i->eps, 0 ); CHKERRQ(ierr);
  ierr = PetscOptionsGetInt( 0, "-iim_Np", &i->Np, 0 ); CHKERRQ(ierr);

  ierr = LocalCoorCreate( &i->lc ); CHKERRQ(ierr);
  ierr = LeastSqCreate(i->Np, is2D, &i->lsq); CHKERRQ(ierr);
  ierr = SpatialIndexCreate( "irregNodes", &i->sidx); CHKERRQ(ierr);

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
  ierr = SpatialIndexDestroy(iim->sidx); CHKERRQ(ierr);
  ierr = ArrayDestroy(iim->coor); CHKERRQ(ierr);
  ierr = ArrayDestroy(iim->val); CHKERRQ(ierr);
  ierr = ArrayDestroy(iim->idx); CHKERRQ(ierr);
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

PetscErrorCode IIMSetEps( IIM iim, PetscReal eps )
{
  iim->eps = eps;
  return 0;
}

PetscErrorCode IIMSetViscosity( IIM iim, PetscReal mu )
{
  iim->mu = mu;
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
