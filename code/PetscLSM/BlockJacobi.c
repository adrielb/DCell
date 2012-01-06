#include "petscksp.h"

PetscErrorCode Seq(PetscInt m);
PetscErrorCode Mpi(PetscInt m);

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc, char **args)
{
  PetscErrorCode  ierr;
  PetscInitialize(&argc, &args, (char *)0, "");
  PetscPrintf(PETSC_COMM_WORLD, "Start\n");
  int m = 700;
  PetscOptionsGetInt(PETSC_NULL, "-m", &m, PETSC_NULL);
  Mpi(m);
  PetscPrintf(PETSC_COMM_WORLD, "End\n");
  ierr = PetscFinalize(); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "Mpi"
PetscErrorCode Mpi(PetscInt m)
{
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  PetscPrintf(PETSC_COMM_WORLD,"Mpi(%d)\n",m);
  Mat M;
  MatCreate(PETSC_COMM_WORLD, &M);
  MatSetSizes(M, PETSC_DECIDE,PETSC_DECIDE, m*m,m*m);
  MatSetType(M, MATMPISBAIJ);
  MatSetOption(M, MAT_SYMMETRIC);
  MatSetOption(M, MAT_ROWS_SORTED);
//  MatSetOption(M,MAT_IGNORE_LOWER_TRIANGULAR);
  
  PetscInt begin, end;
  MatGetOwnershipRange(M,&begin,&end);
  PetscPrintf(PETSC_COMM_SELF,"(%d, %d)\n", begin, end);
  
  int j;
  const PetscInt  COLS[3]={0, 1, m};
  PetscInt cols[3];
  PetscReal vals[3]={4,-1,-1};
  for( int i = begin; i < (end==m*m?end-m:end); ++i)
  {
    for( j=0; j < 3; j++ )
      cols[j] = COLS[j] + i;
    MatSetValues(M, 1, &i, 3, cols, vals, INSERT_VALUES );
  }
  for( int i = end-m; i < end-1; i++ )
  {
    for( j=0; j < 2; j++ )
      cols[j] = COLS[j] + i;
    MatSetValues(M, 1, &i, 2, cols, vals, INSERT_VALUES );
  }
  MatSetValue(M,end-1,end-1,vals[0],INSERT_VALUES);
  
  
  MatAssemblyBegin(M, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(M, MAT_FINAL_ASSEMBLY);
//  ierr = PetscViewerSetFormat(PETSC_VIEWER_STDOUT_SELF, PETSC_VIEWER_ASCII_DENSE); CHKERRQ(ierr);
//  MatView(M,0);
  
  Vec x, b;
  VecCreate(PETSC_COMM_WORLD, &b);
  VecSetSizes(b, PETSC_DECIDE,m*m);
  VecSetType(b,VECMPI);
  VecDuplicate(b,&x);
  VecSet(b, 1./(m*m));
  
  KSP ksp;
  PC pc;
  KSPCreate(PETSC_COMM_WORLD,&ksp);
  KSPSetType(ksp,KSPCG);
  KSPGetPC(ksp, &pc);
  PCSetType(pc,PCBJACOBI);
  KSPSetOperators(ksp,M,M,DIFFERENT_NONZERO_PATTERN);
  KSPSetMonitor(ksp,KSPDefaultMonitor,PETSC_NULL,PETSC_NULL);
  
  KSPSetUp(ksp);
  KSP *subksp;
  PC subpc;
  PCBJacobiGetSubKSP(pc,PETSC_NULL,PETSC_NULL, &subksp);
  KSPSetType(subksp[0],KSPPREONLY);
  KSPGetPC(subksp[0],&subpc);
  PCSetType(subpc,PCICC);
  PCFactorSetLevels(subpc,4);
  Mat mat;
  KSPGetOperators(subksp[0],&mat,PETSC_NULL,PETSC_NULL);
//  PetscPrintf(PETSC_COMM_SELF,"-------------============");
//  MatView(mat,PETSC_VIEWER_STDOUT_SELF);
  
  KSPSolve(ksp,b,x);
//  VecView(x,PETSC_VIEWER_STDOUT_WORLD);
  KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD);
  
  KSPDestroy(ksp);
  VecDestroy(b);
  VecDestroy(x);
  MatDestroy(M);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "Seq"
PetscErrorCode Seq(PetscInt m)
{
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  PetscPrintf(PETSC_COMM_WORLD,"Seq(%d)\n",m);
  Mat M;
  MatCreate(PETSC_COMM_SELF, &M);
  MatSetSizes(M, PETSC_DECIDE,PETSC_DECIDE, m*m,m*m);
  ierr = MatSetType(M, MATSEQSBAIJ); CHKERRQ(ierr);
  
  MatSetOption(M, MAT_SYMMETRIC);
  MatSetOption(M, MAT_ROWS_SORTED);
//  MatSetOption(M,MAT_IGNORE_LOWER_TRIANGULAR);
  
  int j;
  const PetscInt  COLS[3]={0, 1, m};
  PetscInt cols[3];
  PetscReal vals[3]={4,-1,-1};
  for( int i = 0; i < m*m-m; ++i)
  {
    for( j=0; j < 3; j++ )
      cols[j] = COLS[j] + i;
    ierr = MatSetValues(M, 1, &i, 3, cols, vals, INSERT_VALUES ); CHKERRQ(ierr);
  }
  for( int i = m*m-m; i < m*m-1; i++ )
  {
    for( j=0; j < 2; j++ )
      cols[j] = COLS[j] + i;
    ierr = MatSetValues(M, 1, &i, 2, cols, vals, INSERT_VALUES ); CHKERRQ(ierr);
  }
  ierr = MatSetValue(M,m*m-1,m*m-1,vals[0],INSERT_VALUES); CHKERRQ(ierr);
  
  MatAssemblyBegin(M, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(M, MAT_FINAL_ASSEMBLY);
//  ierr = PetscViewerSetFormat(PETSC_VIEWER_STDOUT_SELF, PETSC_VIEWER_ASCII_DENSE); CHKERRQ(ierr);
//  MatView(M,0);
  
  Vec x, b;
  VecCreate(PETSC_COMM_SELF, &b);
  VecSetSizes(b, PETSC_DECIDE,m*m);
  VecSetType(b,VECSEQ);
  VecDuplicate(b,&x);
  VecSet(b, 1./(m*m));
  
  KSP ksp;
  PC pc;
  KSPCreate(PETSC_COMM_WORLD,&ksp);
  KSPSetType(ksp,KSPCG);
  KSPGetPC(ksp, &pc);
  PCSetType(pc,PCICC);
  PCFactorSetLevels(pc,4);
  KSPSetOperators(ksp,M,M,DIFFERENT_NONZERO_PATTERN);
  KSPSetMonitor(ksp,KSPDefaultMonitor,PETSC_NULL,PETSC_NULL);
  
  KSPSolve(ksp,b,x);
//  VecView(x,PETSC_VIEWER_STDOUT_WORLD);
  KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD);
  
  KSPDestroy(ksp);
  VecDestroy(b);
  VecDestroy(x);
  MatDestroy(M);
  PetscFunctionReturn(0);
}