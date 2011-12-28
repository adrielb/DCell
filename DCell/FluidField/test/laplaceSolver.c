#include "FluidField.h"

PetscErrorCode KSPGmres( Mat mat, KSP *k );
PetscErrorCode KSPMG( Mat mat, KSP *k );
PetscErrorCode KSPSOR( Mat mat, KSP *k );
PetscErrorCode KSPICC( Mat mat, KSP *k );
PetscErrorCode KSPASM( Mat mat, KSP *k );
PetscErrorCode KSPCholesky( Mat mat, KSP *k );
PetscErrorCode LaplaceAssemble_2D( Array dbc, DA da, Coor dH, Mat *m );
PetscErrorCode IdentityAssemble_2D( Array dbc, DA da, Coor dH, Mat *m );
PetscErrorCode EnforceBC( DA da, Vec b, Array dbC );

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc, char **args)
{
  PetscErrorCode  ierr;
  ierr = PetscInitialize(&argc, &args, (char *) 0, ""); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "Start %s\n", __FILE__); CHKERRQ(ierr);

  int dx = 512;
  ierr = PetscOptionsGetInt("","-dx",&dx,0); CHKERRQ(ierr);
  Coor dh = {1./dx,1./dx,0};
  int dof = 1;
  iCoor dims = {dx,dx,0};
  DA da;
  ierr = DACreate2d(PETSC_COMM_WORLD,DA_NONPERIODIC,DA_STENCIL_STAR,
              dims.x,dims.y, PETSC_DECIDE,PETSC_DECIDE, dof,1, 0,0, &da); CHKERRQ(ierr);
  Array dbc;
  ierr = ArrayCreate("dbc",sizeof(MatStencil),1000,&dbc); CHKERRQ(ierr);
  Mat mat;
  ierr = LaplaceAssemble_2D(dbc,da,dh,&mat); CHKERRQ(ierr);
  KSP ksp;
//  ierr = KSPGmres( mat, &ksp ); CHKERRQ(ierr);
//  ierr = KSPMG( mat, &ksp ); CHKERRQ(ierr);
//  ierr = KSPCholesky( mat, &ksp ); CHKERRQ(ierr);
  ierr = KSPICC( mat, &ksp ); CHKERRQ(ierr);
//  ierr = KSPSOR( mat, &ksp ); CHKERRQ(ierr);
//  ierr = KSPASM( mat, &ksp ); CHKERRQ(ierr);
  Vec x,b;
  ierr = DACreateGlobalVector(da,&x); CHKERRQ(ierr);
  ierr = DACreateGlobalVector(da,&b); CHKERRQ(ierr);
  ierr = VecSet(b,1.); CHKERRQ(ierr);
  ierr = EnforceBC(da,b,dbc); CHKERRQ(ierr);
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
  ierr = ArrayDestroy(dbc); CHKERRQ(ierr);
  ierr = MatDestroy(mat); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "End %s\n", __FILE__); CHKERRQ(ierr);
  ierr = PetscFinalize(); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

// src/ksp/ksp/examples/tests/ex19.c

PetscErrorCode KSPMG( Mat mat, KSP *k )
{
  int i;
  PetscErrorCode ierr;
  KSP ksp;
  ierr = KSPCreate(PETSC_COMM_WORLD,&ksp); CHKERRQ(ierr);
  ierr = KSPSetOperators(ksp,mat,mat,SAME_PRECONDITIONER); CHKERRQ(ierr);
  ierr = KSPSetType(ksp,KSPGMRES); CHKERRQ(ierr);
  ierr = KSPSetFromOptions(ksp); CHKERRQ(ierr);

  PC pc;
  ierr = KSPGetPC(ksp,&pc); CHKERRQ(ierr);
  ierr = PCSetType(pc,PCMG); CHKERRQ(ierr);
  int nlevels = 2;
  MPI_Comm *comms;
  ierr  = PetscMalloc(nlevels*sizeof(MPI_Comm),&comms);CHKERRQ(ierr);
  for ( i=0; i<nlevels; i++) {
    comms[i] = PETSC_COMM_WORLD;
  }
  ierr = PCMGSetLevels(pc,nlevels,comms); CHKERRQ(ierr);
  ierr = PetscFree(comms);CHKERRQ(ierr);
  ierr = PCMGSetType(pc, PC_MG_FULL); CHKERRQ(ierr);

  DA dac;
  iCoor dimsf = {128,128,0};
  iCoor dims = { ceil((dimsf.x+1)/2), ceil((dimsf.y+1)/2), 0 };
  int dof = 1;
  ierr = DACreate2d(PETSC_COMM_WORLD,DA_NONPERIODIC,DA_STENCIL_STAR,
              dims.x,dims.y, PETSC_DECIDE,PETSC_DECIDE, dof,1, 0,0, &dac); CHKERRQ(ierr);

  DA *das;
  ierr = PetscMalloc(nlevels*sizeof(DA),&das);CHKERRQ(ierr);
  ierr = DARefineHierarchy(dac,nlevels,das);CHKERRQ(ierr);
  Mat R;
  for ( i=1; i<nlevels; i++) {
    ierr = DAGetInterpolation(das[i-1],das[i],&R,0); CHKERRQ(ierr);
    ierr = PCMGSetInterpolation(pc,i,R);CHKERRQ(ierr);
    ierr = PCMGSetRestriction(pc,i,R);CHKERRQ(ierr);
  }

  /* Create coarse level */
//  ierr = PCMGGetCoarseSolve(pc,&user.ksp_coarse);CHKERRQ(ierr);
//  ierr = KSPSetOptionsPrefix(user.ksp_coarse,"coarse_");CHKERRQ(ierr);
//  ierr = KSPSetFromOptions(user.ksp_coarse);CHKERRQ(ierr);
//  ierr = KSPSetOperators(user.ksp_coarse,user.coarse.J,user.coarse.J,DIFFERENT_NONZERO_PATTERN);CHKERRQ(ierr);
//  ierr = PCMGSetX(pc,COARSE_LEVEL,user.coarse.x);CHKERRQ(ierr);
//  ierr = PCMGSetRhs(pc,COARSE_LEVEL,user.coarse.b);CHKERRQ(ierr);

  ierr = PCSetUp(pc); CHKERRQ(ierr);

  *k = ksp;
  return 0;
}

PetscErrorCode KSPASM( Mat mat, KSP *k )
{
  PetscErrorCode ierr;
  KSP ksp;
  ierr = KSPCreate(PETSC_COMM_WORLD,&ksp); CHKERRQ(ierr);
  ierr = KSPSetOperators(ksp,mat,mat,SAME_PRECONDITIONER); CHKERRQ(ierr);
  ierr = KSPSetType(ksp,KSPGMRES); CHKERRQ(ierr);
  ierr = KSPSetTolerances(ksp,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT); CHKERRQ(ierr);
  ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);

  PC pc;
  ierr = KSPGetPC(ksp,&pc); CHKERRQ(ierr);
  ierr = PCSetType(pc,PCASM); CHKERRQ(ierr);
  ierr = PCASMSetOverlap(pc,50); CHKERRQ(ierr);
  ierr = PCSetUp(pc); CHKERRQ(ierr);
  KSP *subksp;
  ierr = PCASMGetSubKSP(pc,PETSC_NULL,PETSC_NULL,&subksp); CHKERRQ(ierr);
  ierr = KSPSetType(subksp[0],KSPPREONLY); CHKERRQ(ierr);
  ierr = KSPGetPC(subksp[0],&pc); CHKERRQ(ierr);
  ierr = PCSetType(pc,PCCHOLESKY); CHKERRQ(ierr);
  ierr = PCFactorSetMatOrderingType(pc,MATORDERING_ND); CHKERRQ(ierr);
  ierr = PCSetUp(pc); CHKERRQ(ierr);
  *k = ksp;
  return 0;
}

PetscErrorCode KSPSOR( Mat mat, KSP *k )
{
  PetscErrorCode ierr;
  KSP ksp;
  ierr = KSPCreate(PETSC_COMM_WORLD,&ksp); CHKERRQ(ierr);
  ierr = KSPSetOperators(ksp,mat,mat,SAME_PRECONDITIONER); CHKERRQ(ierr);
  ierr = KSPSetType(ksp,KSPGMRES); CHKERRQ(ierr);
  ierr = KSPSetTolerances(ksp,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT); CHKERRQ(ierr);

  PC pc;
  ierr = KSPGetPC(ksp,&pc); CHKERRQ(ierr);
  ierr = PCSetType(pc,PCSOR); CHKERRQ(ierr);
  ierr = PCSORSetIterations(pc,1,1); CHKERRQ(ierr);
  ierr = PCSORSetSymmetric(pc,SOR_SYMMETRIC_SWEEP); CHKERRQ(ierr);
  ierr = PCSORSetOmega(pc,1.9); CHKERRQ(ierr);
  ierr = PCSetUp(pc); CHKERRQ(ierr);
  *k = ksp;
  return 0;
}

PetscErrorCode KSPICC( Mat mat, KSP *k )
{
  PetscErrorCode ierr;
  KSP ksp;
  ierr = KSPCreate(PETSC_COMM_WORLD,&ksp); CHKERRQ(ierr);
  ierr = KSPSetOperators(ksp,mat,mat,SAME_PRECONDITIONER); CHKERRQ(ierr);
  ierr = KSPSetType(ksp,KSPGMRES); CHKERRQ(ierr);
//  ierr = KSPSetInitialGuessNonzero(ksp,PETSC_TRUE); CHKERRQ(ierr);
  PC pc;
  ierr = KSPGetPC(ksp,&pc); CHKERRQ(ierr);
  ierr = PCSetType(pc,PCICC); CHKERRQ(ierr);
  ierr = PCFactorSetMatOrderingType(pc,MATORDERING_ND); CHKERRQ(ierr);
  ierr = PCFactorSetLevels(pc, 6000); CHKERRQ(ierr);
  ierr = PCFactorSetFill(pc,5.805); CHKERRQ(ierr);
  ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);
  PetscLogDouble t1,t2;
  ierr = PetscGetTime(&t1); CHKERRQ(ierr);
  ierr = PCSetUp(pc); CHKERRQ(ierr);
  ierr = PetscGetTime(&t2); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"%f sec\n",t2-t1); CHKERRQ(ierr);
  *k = ksp;
  return 0;
}

PetscErrorCode KSPCholesky( Mat mat, KSP *k )
{
  PetscErrorCode ierr;
  KSP ksp;
  ierr = KSPCreate(PETSC_COMM_WORLD,&ksp); CHKERRQ(ierr);
  ierr = KSPSetOperators(ksp,mat,mat,SAME_PRECONDITIONER); CHKERRQ(ierr);
  ierr = KSPSetType(ksp,KSPPREONLY); CHKERRQ(ierr);
  PC pc;
  ierr = KSPGetPC(ksp,&pc); CHKERRQ(ierr);
  ierr = PCSetType(pc,PCCHOLESKY); CHKERRQ(ierr);
  ierr = PCFactorSetMatOrderingType(pc,MATORDERING_ND); CHKERRQ(ierr);
//  ierr = PCFactorSetFill(pc,5.805); CHKERRQ(ierr);
  ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);
  PetscLogDouble t1,t2;
  ierr = PetscGetTime(&t1); CHKERRQ(ierr);
  ierr = PCSetUp(pc); CHKERRQ(ierr);
  ierr = PetscGetTime(&t2); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"%f sec\n",t2-t1); CHKERRQ(ierr);
  *k = ksp;
  return 0;
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
#define __FUNCT__ "LaplaceAssemble_2D"
PetscErrorCode LaplaceAssemble_2D( Array dbc, DA da, Coor dH, Mat *m )
{
  DALocalInfo i;
  Mat mat;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  MPI_Comm comm;
  PetscObjectGetComm( (PetscObject)da, &comm);
  int size;
  MPI_Comm_size(comm, &size);
  if( size == 1 ) {
    ierr = DAGetMatrix(da, MATSEQAIJ, &mat); CHKERRQ(ierr);
  } else {
    ierr = DAGetMatrix(da, MATMPIAIJ, &mat); CHKERRQ(ierr);
  }

  DAGetLocalInfo(da, &i);
  PetscReal hx = dH.x, hy = dH.y;
  const int n = 5; // 5-point laplace stencil + 2-point gradient stencil
  MatStencil row, col[5];
  PetscReal val[5] = {-hx*hx,-hx*hx,-hy*hy,-hy*hy,2*(hx*hx+hy*hy)},
            val_D[5] = { 0, 0, 0, 0, 1};
  int c;
  for( row.j = i.ys; row.j < i.ys+i.ym; ++row.j)
  {
    for( row.i = i.xs; row.i < i.xs+i.xm; ++row.i)
    {
      row.c = 0;
      for( c = 0; c < n; ++c)
        PetscMemcpy(&col[c],&row,sizeof(MatStencil));
      col[0].i--; col[1].i++;  // [-hx, hx]
      col[2].j--; col[3].j++;  // [-hy, hy]
      if( col[1].i >= i.mx ) col[1].i = -1;
      if( col[3].j >= i.my ) col[3].j = -1;

      // Laplace(u)
      if( 0 < row.i && row.i < i.mx - 1 &&
          0 < row.j && row.j < i.my - 2 ) {
        ierr = MatSetValuesStencil(mat,1,&row,n,col,val,INSERT_VALUES); CHKERRQ(ierr);
      } else {
        ierr = MatSetValuesStencil(mat,1,&row,5,col,val_D,INSERT_VALUES); CHKERRQ(ierr);
        ierr = FluidField_AppendDBC(dbc,row); CHKERRQ(ierr);
      }
    } // row.i
    ierr = PetscPrintf(PETSC_COMM_WORLD,"."); CHKERRQ(ierr);
  } // row.j
  ierr = PetscPrintf(PETSC_COMM_WORLD,"\n"); CHKERRQ(ierr);

  ierr = MatAssemblyBegin(mat,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(mat,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
  /*
   *  Make velocity block symmetric by eliminating BCs
   */
  for( row.j = i.ys; row.j < i.ys+i.ym; ++row.j)
  {
    for( row.i = i.xs; row.i < i.xs+i.xm; ++row.i)
    {
      row.c = 0;
      for( c = 0; c < 5; ++c)
        PetscMemcpy(&col[c],&row,sizeof(MatStencil));
      col[0].i--; col[1].i++;  // [-hx, hx]
      col[2].j--; col[3].j++;  // [-hy, hy]
      if( col[1].i >= i.mx ) col[1].i = -1;
      if( col[3].j >= i.my ) col[3].j = -1;

      // Laplace(u)
      if( 0 < row.i && row.i < i.mx - 1 &&
          0 < row.j && row.j < i.my - 2 ) {

      } else {
        ierr = MatSetValuesStencil(mat,5,col,1,&row,val_D,INSERT_VALUES); CHKERRQ(ierr);
      }
    } // row.i
    ierr = PetscPrintf(PETSC_COMM_WORLD,"."); CHKERRQ(ierr);
  } // row.j
  ierr = PetscPrintf(PETSC_COMM_WORLD,"\n"); CHKERRQ(ierr);

  ierr = MatAssemblyBegin(mat,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(mat,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

  *m = mat;

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "IdentityAssemble_2D"
PetscErrorCode IdentityAssemble_2D( Array dbc, DA da, Coor dH, Mat *m )
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

#undef __FUNCT__
#define __FUNCT__ "EnforceBC"
PetscErrorCode EnforceBC( DA da, Vec b, Array dbC )
{
  int i;
  int len = ArrayLength(dbC);
  MatStencil *dbc;
  PetscReal ***rhs;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = DAVecGetArrayDOF(da,b,&rhs); CHKERRQ(ierr);
  dbc = ArrayGetData(dbC);
  for ( i = 0; i < len; ++i) {
    rhs[dbc[i].j][dbc[i].i][dbc[i].c] = 0;
  }
  ierr = DAVecRestoreArrayDOF(da,b,&rhs); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
