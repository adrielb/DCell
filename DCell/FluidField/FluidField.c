#include "FluidField.h"
#include "FluidField_private.h"

#undef __FUNCT__
#define __FUNCT__ "FluidFieldCreate"
PetscErrorCode FluidFieldCreate(MPI_Comm comm, FluidField *fluid)
{
  FluidField f;
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  ierr = PetscNew(struct _FluidField, &f); CHKERRQ(ierr);

  // Default values

  f->comm = comm;
  f->mu = 1e-3; // pN s / um^2
  ierr = PetscOptionsGetReal(0,"-fluid_viscosity",&f->mu,0); CHKERRQ(ierr);

  PetscReal dx = 1;
  ierr = PetscOptionsGetReal(0,"-fluid_dx",&dx,0); CHKERRQ(ierr);
  f->dh.x = dx;
  f->dh.y = dx;
  f->dh.z = dx;
  
  f->dims.x = 7;
  f->dims.y = 7;
  f->dims.z = 0;

  f->is3D = PETSC_FALSE;

  // Create BC index set
  ierr = ArrayCreate("dirichletBC",sizeof(MatStencil),128,&f->dirichletBC); CHKERRQ(ierr);

  ierr = FluidFieldRegisterEvents(); CHKERRQ(ierr);

  *fluid = f;

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "FluidFieldDestroy"
PetscErrorCode FluidFieldDestroy(FluidField f)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = DMDestroy(&f->daV); CHKERRQ(ierr);
  ierr = MatDestroy(&f->mat); CHKERRQ(ierr);
  ierr = KSPDestroy(&f->ksp); CHKERRQ(ierr);
  ierr = VecDestroy(&f->rhs); CHKERRQ(ierr);
  ierr = VecDestroy(&f->vel); CHKERRQ(ierr);
  ierr = ArrayDestroy(f->dirichletBC); CHKERRQ(ierr);
  GA_Destroy(f->ga);

  ierr = DMDestroy(&f->daE); CHKERRQ(ierr);
//  ierr = VecDestroy(f->E); CHKERRQ(ierr);

  ierr = DMDestroy(&f->daB); CHKERRQ(ierr);
  ierr = VecDestroy(&f->buf); CHKERRQ(ierr);
  ierr = PetscFree(f); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode FluidFieldSetDims( FluidField f, iCoor dims )
{
  f->dims = dims;
  f->is3D = dims.z > 1;
  return 0;
}

PetscErrorCode FluidFieldSetDx( FluidField f, PetscReal dx )
{
  f->dh.x = dx;
  f->dh.y = dx;
  f->dh.z = f->is3D ? dx : 0;
  return 0;
}

PetscErrorCode FluidFieldSetMask( FluidField f, Grid mask )
{
  f->mask = mask;
  return 0;
}

PetscErrorCode FluidFieldSetViscosity( FluidField f, PetscReal mu )
{
  f->mu = mu;
  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "FluidFieldSetup"
PetscErrorCode FluidFieldSetup( FluidField f )
{
  PetscLogDouble t1,t2;
  PetscErrorCode ierr;

  PetscFunctionBegin;

  // Assemble viscous matricies
  ierr = FluidFieldMatAssemble( f ); CHKERRQ(ierr);

  ierr = PetscGetTime(&t1); CHKERRQ(ierr);

  // Create vectors
  ierr = GACreate( f->daV, &f->ga); CHKERRQ(ierr);
  ierr = DMCreateGlobalVector(f->daV,&f->rhs); CHKERRQ(ierr);
  ierr = VecDuplicate(f->rhs,&f->vel); CHKERRQ(ierr);
//  ierr = DACreateGlobalVector(f->daE,&f->E); CHKERRQ(ierr);
  ierr = DMCreateGlobalVector(f->daB,&f->buf); CHKERRQ(ierr);

  // Set up the outer solver
  ierr = KSPCreate(f->comm,&f->ksp); CHKERRQ(ierr);
  ierr = KSPSetOperators(f->ksp,f->mat,f->mat, SAME_PRECONDITIONER); CHKERRQ(ierr);
  ierr = KSPSetType(f->ksp,KSPFGMRES); CHKERRQ(ierr);
  ierr = KSPSetInitialGuessNonzero(f->ksp,PETSC_TRUE); CHKERRQ(ierr);
  ierr = KSPSetTolerances(f->ksp,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT); CHKERRQ(ierr);
  ierr = KSPSetFromOptions(f->ksp);CHKERRQ(ierr);

  // Split pressure from velocity [ u v w | p ]
  PC pc;
  ierr = KSPGetPC(f->ksp,&pc); CHKERRQ(ierr);
  ierr = PCSetType(pc, PCFIELDSPLIT); CHKERRQ(ierr);
  ierr = PCFieldSplitSetType(pc,PC_COMPOSITE_SCHUR); CHKERRQ(ierr);
  if( f->is3D ) {
    ierr = PCFieldSplitSetBlockSize(pc,4); CHKERRQ(ierr);                            // [p u v w]
    ierr = PCFieldSplitSetFields(pc,"v",3,(int[]){U_FACE,V_FACE,W_FACE}); CHKERRQ(ierr); // [u v w]
    ierr = PCFieldSplitSetFields(pc,"p",1,(int[]){CELL_CENTER}); CHKERRQ(ierr);          // [ p ]
  } else {
    ierr = PCFieldSplitSetBlockSize(pc,3); CHKERRQ(ierr);                      // [p u v]
    ierr = PCFieldSplitSetFields(pc,"v",2,(int[]){U_FACE,V_FACE}); CHKERRQ(ierr);  // [u v]
    ierr = PCFieldSplitSetFields(pc,"p",1,(int[]){CELL_CENTER}); CHKERRQ(ierr);    // [ p ]
  }
  ierr = PCSetUp(pc); CHKERRQ(ierr);
  int nVelP;
  KSP *kspVelP;
  ierr = PCFieldSplitGetSubKSP(pc,&nVelP,&kspVelP); CHKERRQ(ierr);
  ierr = KSPSetTolerances(kspVelP[1],PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,4); CHKERRQ(ierr);
  ierr = KSPSetType(kspVelP[1],KSPGMRES); CHKERRQ(ierr);
  ierr = KSPGetPC(kspVelP[1],&pc); CHKERRQ(ierr);
  ierr = PCSetType(pc,PCNONE); CHKERRQ(ierr);
  ierr = KSPSetFromOptions(kspVelP[1]);CHKERRQ(ierr);

  // Split velocity [u v w] into component matricies [u], [v], [w]
  ierr = KSPSetType(kspVelP[0],KSPPREONLY); CHKERRQ(ierr);
  ierr = KSPGetPC(kspVelP[0],&pc); CHKERRQ(ierr);
  ierr = PCSetType(pc, PCFIELDSPLIT); CHKERRQ(ierr);
  ierr = PCFieldSplitSetType(pc,PC_COMPOSITE_ADDITIVE); CHKERRQ(ierr);
  ierr = PCFieldSplitSetBlockSize(pc,f->is3D?3:2); CHKERRQ(ierr);
  ierr = PCSetUp(pc); CHKERRQ(ierr);

  /* Set solver for each velocity component
   * Split component velocity as parallel blocks along processors
   * Use direct solver for each block
   * TODO: use MG, w/FFT on coarse grid
   */

  ierr = PetscGetTime(&t2); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Finished Solver Setup: %f sec\n",t2-t1); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "FluidFieldWrite"
PetscErrorCode FluidFieldWrite(FluidField f, int t)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = VecWrite(f->vel, "uvp", t); CHKERRQ(ierr);
//  ierr = VecWrite(f->E,   "E",   t); CHKERRQ(ierr); //TODO: uncomment when viscoelastic effects coded
  ierr = VecWrite(f->rhs, "rhs", t); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "FluidFieldRegisterEvents"
static int EVENTS_registered = PETSC_FALSE;
PetscErrorCode FluidFieldRegisterEvents(  )
{
  PetscErrorCode ierr;
  if( EVENTS_registered )
    PetscFunctionReturn(0);

  ierr = PetscLogEventRegister("FluidMaxVel",0,&EVENT_FluidFieldMaxVelocityMag); CHKERRQ(ierr);
  ierr = PetscLogEventRegister("FluidDiscreteCC",0,&EVENT_FluidField_DiscreteCompatibilityCondition); CHKERRQ(ierr);


  EVENTS_registered = PETSC_TRUE;
  PetscFunctionReturn(0);
}
