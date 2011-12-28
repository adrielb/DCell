#include "petsc.h"

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc, char **args)
{
  PetscErrorCode  ierr;
  ierr = PetscInitialize(&argc, &args, (char *) 0, ""); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "Start %s\n", __FILE__); CHKERRQ(ierr);
  MPI_Comm comm = PETSC_COMM_WORLD;
  int n = 10;
  Mat mat;
  ierr = MatCreate(comm, &mat);
  ierr = MatSetType(mat, MATSEQAIJ);
  ierr = MatSetSizes(mat, PETSC_DECIDE,PETSC_DECIDE, n, n);
  int i;
  for( i = 0; i < n; i++) {
    int col[] = {i};
    int row[] = {i};
    PetscReal val = n;
    MatSetValues( mat, 1, col, 1, row, &val, INSERT_VALUES);
  }
  ierr = MatAssemblyBegin(mat, MAT_FINAL_ASSEMBLY);
  ierr = MatAssemblyEnd(  mat, MAT_FINAL_ASSEMBLY);
  
  KSP ksp;
  ierr = KSPCreate( comm, &ksp );
  ierr = KSPSetType( ksp, KSPGMRES );
  ierr = KSPSetFromOptions(ksp);
  PC pc;
  ierr = KSPGetPC( ksp, &pc);
  ierr = PCSetType( pc, PCNONE);
  ierr = KSPSetOperators( ksp, mat, mat, DIFFERENT_NONZERO_PATTERN);
  Vec b, x;
  ierr = VecCreateSeq( comm, n, &b);
  ierr = VecCreateSeq( comm, n, &x);
  ierr = VecSet(b,1);
  ierr = KSPSolve(ksp, b, x);
  ierr = KSPDestroy(ksp);
  ierr = VecDestroy(b);
  ierr = VecDestroy(x);
  ierr = MatDestroy(mat); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "End %s\n", __FILE__); CHKERRQ(ierr);
  ierr = PetscFinalize(); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
