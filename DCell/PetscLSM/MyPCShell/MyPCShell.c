#include "petscksp.h"
#include "petscdmmg.h"
#include "petscmg.h"

PetscErrorCode MyPCShellApply( void* ctx, Vec, Vec );

typedef struct {
	
} MyPCShellContext;
PetscInt EVENT_MyPCShellApply;

PetscErrorCode ComputeRHS(DMMG dmmg,Vec b) {return 0;}
PetscErrorCode ComputeJacobian(DMMG dmmg,Mat jac,Mat B)
{
  /*
  Vec vOne;
  VecCreate( PETSC_COMM_SELF, vOne );
  VecSetSizes(vOne, len, len);
  VecSetType(vOne, VECSEQ);
  VecSet(vOne, 1.);
  MatDiagonalSet(B, vOne, INSERT_VALUES);
  */
  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc, char **argv)
{
	PetscErrorCode ierr;
	KSP ksp;
	PC pc, pcmg;
	Vec b, x;
	Mat A;
	DMMG *dmmg;
  DA da;
	
	
	PetscInitialize(&argc, &argv, (char*) 0, "");
	PetscLogEventRegister(&EVENT_MyPCShellApply,"MyPCShellApply", 0);
	
	VecCreate(PETSC_COMM_SELF, &b);
	VecSetSizes(b, 10, 10);
	VecSetType(b, VECSEQ);
	VecDuplicate(b, &x);
	
	MatCreate(PETSC_COMM_SELF, &A);
	MatSetSizes(A,10,10,10,10);
	MatSetType(A, MATSEQAIJ);
  
  DMMGCreate(PETSC_COMM_SELF, 3, PETSC_NULL, &dmmg);
  DACreate1d(PETSC_COMM_SELF, DA_XPERIODIC, 5, 1, 1, PETSC_NULL, &da);
  DMMGSetDM( dmmg, (DM) da );
  DADestroy(da);
  
  ierr = DMMGSetKSP(dmmg, ComputeRHS, ComputeJacobian); CHKERRQ(ierr);
  ierr = DMMGSetUp(dmmg); CHKERRQ(ierr);

  ksp = DMMGGetKSP(dmmg);
  KSPSetType(ksp, KSPPREONLY);
  KSPGetPC(ksp, &pcmg);
  for( int i = 0; i < DMMGGetLevels(dmmg); i++)
  {
  	ierr = PCMGGetSmoother(pcmg, i, &ksp); CHKERRQ(ierr);
  	KSPSetType(ksp, KSPRICHARDSON);
		KSPGetPC(ksp, &pc);
		PCSetType(pc, PCSHELL);
		PCShellSetApply(pc, MyPCShellApply);
  }
  
  ierr = DMMGView(dmmg, PETSC_VIEWER_STDOUT_SELF); CHKERRQ(ierr);
  ierr = DMMGSolve(dmmg); CHKERRQ(ierr);
  
/*
  KSPCreate(PETSC_COMM_SELF, &ksp);

	
	KSPSetOperators(ksp, A, A, SAME_NONZERO_PATTERN);
	KSPSetMonitor(ksp, KSPDefaultMonitor, PETSC_NULL, PETSC_NULL);
	ierr = KSPSolve(ksp, b, x); CHKERRQ(ierr);
*/	
	ksp = DMMGGetKSP(dmmg);
	PetscInt iter;
	KSPGetIterationNumber(ksp, &iter);
	PetscPrintf(PETSC_COMM_SELF, "numiter: %d\n",iter);
	KSPConvergedReason reason;
	KSPGetConvergedReason(ksp, &reason);
	PetscPrintf(PETSC_COMM_SELF, "reason: %d\n", reason);

	VecDestroy(b);
	VecDestroy(x);
	DMMGDestroy(dmmg);
	PetscFinalize();
}

#undef __FUNCT__
#define __FUNCT__ "MyPCShellApply"
PetscErrorCode MyPCShellApply( void* ctx, Vec x, Vec y)
{
	PetscErrorCode ierr;
	PetscInt size;
	
	PetscFunctionBegin;
	PetscLogEventBegin(EVENT_MyPCShellApply,0,0,0,0);
  
	VecGetSize(x, &size);
  PetscPrintf(PETSC_COMM_WORLD, "size: %d\n", size );
//	VecView(x, PETSC_VIEWER_STDOUT_SELF);
	
	PetscLogEventEnd(EVENT_MyPCShellApply,0,0,0,0);
	PetscFunctionReturn(0);
}