#include "Common.h"

PetscErrorCode KSPGmres( Mat mat, KSP *k );

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc, char **args)
{
  PetscErrorCode  ierr;
  ierr = PetscInitialize(&argc, &args, (char *) 0, ""); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "Start %s\n", __FILE__); CHKERRQ(ierr);

  int dx = 64;
  ierr = PetscOptionsGetInt("","-dx",&dx,0); CHKERRQ(ierr);
  int dof = 1;
  DA da;
  ierr = DACreate2d(PETSC_COMM_WORLD,DA_NONPERIODIC,DA_STENCIL_STAR,
              dx,dx, PETSC_DECIDE,PETSC_DECIDE, dof,1, 0,0, &da); CHKERRQ(ierr);

  Mat mat;
  IdentityAssemble_2D( da, &mat);
  KSP ksp;
  ierr = KSPGmres( mat, &ksp ); CHKERRQ(ierr);

  Vec x,b;
  ierr = DACreateGlobalVector(da,&x); CHKERRQ(ierr);
  ierr = DACreateGlobalVector(da,&b); CHKERRQ(ierr);
//  ierr = VecCreateSeq(PETSC_COMM_WORLD, dx*dx,&x); CHKERRQ(ierr);
//  ierr = VecCreateSeq(PETSC_COMM_WORLD, dx*dx,&b); CHKERRQ(ierr);
  ierr = VecSet(b,1.); CHKERRQ(ierr);

  PetscLogDouble t1,t2,dt,tavg=0;
  int i, N = 20;
  for ( i = 0; i < N; ++i) {
    ierr = PetscGetTime(&t1); CHKERRQ(ierr);
    ierr = KSPSolve(ksp,b,x); CHKERRQ(ierr);
    ierr = PetscGetTime(&t2); CHKERRQ(ierr);
    dt = t2-t1;
    PetscPrintf(PETSC_COMM_WORLD,"dt: %e\n",dt);
    tavg += dt/N;
  }
  ierr = PetscPrintf(PETSC_COMM_WORLD,"AVG: %f sec\n",tavg); CHKERRQ(ierr);

  ierr = KSPDestroy(ksp); CHKERRQ(ierr);
  ierr = VecDestroy(x); CHKERRQ(ierr);
  ierr = VecDestroy(b); CHKERRQ(ierr);
  ierr = DADestroy(da); CHKERRQ(ierr);
  ierr = MatDestroy(mat); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "End %s\n", __FILE__); CHKERRQ(ierr);
  ierr = PetscFinalize(); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


PetscErrorCode KSPGmres( Mat mat, KSP *k )
{
  PetscErrorCode ierr;
  KSP ksp;
  ierr = KSPCreate(PETSC_COMM_WORLD,&ksp); CHKERRQ(ierr);
  ierr = KSPSetOperators(ksp,mat,mat,SAME_PRECONDITIONER); CHKERRQ(ierr);
  ierr = KSPSetType(ksp,KSPGMRES); CHKERRQ(ierr);
  PC pc;
  ierr = KSPGetPC(ksp,&pc); CHKERRQ(ierr);
  ierr = PCSetType(pc,PCNONE); CHKERRQ(ierr);
  *k = ksp;
  return 0;
}


#undef __FUNCT__
#define __FUNCT__ "IdentityAssemble_2D"
PetscErrorCode IdentityAssemble_2D( DA da, Mat *m )
{
  DALocalInfo i;
  Mat mat;
  MatStencil row;
  PetscReal val_D[] = { 1 };
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = DAGetMatrix(da, MATSEQAIJ, &mat); CHKERRQ(ierr);
  ierr = DAGetLocalInfo(da, &i); CHKERRQ(ierr);
  for( row.j = i.ys; row.j < i.ys+i.ym; ++row.j)
  {
    for( row.i = i.xs; row.i < i.xs+i.xm; ++row.i)
    {
      row.c = 0;
      ierr = MatSetValuesStencil(mat,1,&row,1,&row,val_D,INSERT_VALUES); CHKERRQ(ierr);
    } // row.i
    ierr = PetscPrintf(PETSC_COMM_WORLD,"0"); CHKERRQ(ierr);
  } // row.j
  ierr = PetscPrintf(PETSC_COMM_WORLD,"\n"); CHKERRQ(ierr);

  ierr = MatAssemblyBegin(mat,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(mat,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

  *m = mat;

  PetscFunctionReturn(0);
}

void func(DA da)
{
  int ga;

//  GACreate(da, &ga);
  ga = GA_Create_handle();
  GA_Destroy(ga);
}

