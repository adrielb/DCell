#include "MicrofluidicSimulator.h"

#undef __FUNCT__
#define __FUNCT__ "DiffusionAdvection"
PetscInt EVENT_DiffusionAdvection;
PetscErrorCode DiffusionAdvection(  )
{
  PetscErrorCode ierr;
  
  
  PetscFunctionBegin;
  PetscLogEventBegin(EVENT_DiffusionAdvection,0,0,0,0);
  PetscLogEventRegister(&EVENT_DiffusionAdvection,"DiffusionAdvection", 0);
  
   
  /*  SS-Conc Profile  */
  ierr = AssembleSSConcentration(uc); CHKERRQ(ierr);
  ierr = AssembleConcentrationRHS( uc ); CHKERRQ(ierr);
  ierr = KSPSetType( uc->ksp, KSPPREONLY); CHKERRQ(ierr);
  PC pc;
  ierr = KSPGetPC(uc->ksp, &pc); CHKERRQ(ierr);
  ierr = PCSetType(pc, PCLU); CHKERRQ(ierr);
  Mat umf;
  MatConvert(uc->A, MATUMFPACK, MAT_INITIAL_MATRIX,&umf);
  KSPSetOperators(uc->ksp, umf, umf, DIFFERENT_NONZERO_PATTERN);
  ierr = KSPSolve(uc->ksp, uc->b, uc->c); CHKERRQ(ierr);
  ierr = WriteResult(uc, uc->c, "c_concentration"); CHKERRQ(ierr);
  MatDestroy(umf);
PetscFunctionReturn(0); 
  
  /*  Assemble and Solve for Diffusion-Advection  */
ierr = PetscPrintf(PETSC_COMM_WORLD, "***********************\nADVECTION-DIFFUSION\n"); CHKERRQ(ierr);
  ierr = SolveConcentration( uc ); CHKERRQ(ierr);
  ierr = AssembleDiffusionMatrix( uc ); CHKERRQ(ierr);
  ierr = AssembleConcentrationRHS( uc ); CHKERRQ(ierr);
//MatView(uc->A, PETSC_VIEWER_STDOUT_SELF);
  ierr = KSPSetType( uc->ksp, KSPGMRES); CHKERRQ(ierr);
  ierr = KSPSetInitialGuessNonzero(uc->ksp, PETSC_TRUE); CHKERRQ(ierr);
  WriteVec wv;  char filename[30];
  ierr = WriteVecCreate(uc->c, "conc", &wv); CHKERRQ(ierr);
  sprintf( filename,"conc.%d", 0);
  ierr = WriteResult(uc, uc->b, filename); CHKERRQ(ierr);
  ierr = VecCopy(uc->b,uc->c); CHKERRQ(ierr);
  int count = 0;
  for( int i = 1; i < 1000; i++)
  {
    ierr = KSPSolve(uc->ksp, uc->b, uc->c); CHKERRQ(ierr);
//    ierr = WriteVecToDisk(wv); CHKERRQ(ierr);
    if( i % 100 == 0 )
    {
      ierr = PetscPrintf(PETSC_COMM_WORLD, "C: %d\n", count); CHKERRQ(ierr);
      sprintf( filename,"conc.%d", count);
      ierr = WriteResult(uc, uc->c, filename); CHKERRQ(ierr);
      count++;
    }
    ierr = VecSwap(uc->b, uc->c); CHKERRQ(ierr);
  }
  ierr = WriteVecDestroy(wv); CHKERRQ(ierr);
  
  PetscLogEventEnd(EVENT_DiffusionAdvection,0,0,0,0);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "AssembleDiffusionMatrix"
PetscInt EVENT_AssembleDiffusionMatrix;
PetscErrorCode AssembleDiffusionMatrix( UserContext *uc )
{
  int i, j;
  PetscReal *u, *v, ut, vt, val[5];
  Node *n;
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  PetscLogEventBegin(EVENT_AssembleDiffusionMatrix,0,0,0,0);
  VecGetArray(uc->u, &u);
  VecGetArray(uc->v, &v);
  
  for (i = 0; i < uc->numNodes; ++i)
  {
    n = &uc->nodes[i];

    // Diffusive part
    for( j = 0; j < 4; j++)
      val[j] = (-1 * DIFFUSION * DELTA_T) / ( DELTA_X * DELTA_X );
    val[4] = 1. + ( n->numNei * DIFFUSION * DELTA_T ) / ( DELTA_X * DELTA_X );

    // Advective part (upwinding) 
    ut = u[i] * DELTA_T / DELTA_X;
    vt = v[i] * DELTA_T / DELTA_X;
    
    if( u[i] < 0. ) { val[2] += ut; val[4] -= ut; }
    if( u[i] > 0. ) { val[1] -= ut; val[4] += ut; }
    if( v[i] < 0. ) { val[3] += vt; val[4] -= vt; }
    if( v[i] > 0. ) { val[0] -= vt; val[4] += vt; }

    ierr = MatSetValues(uc->A, 1, &i, 5, n->star, val, INSERT_VALUES); CHKERRQ(ierr);
  }
  
  VecRestoreArray(uc->v, &v);
  VecRestoreArray(uc->u, &u);
  PetscLogEventEnd(EVENT_AssembleDiffusionMatrix,0,0,0,0);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "AssembleConcentrationRHS"
PetscInt EVENT_AssembleConcentrationRHS;
PetscErrorCode AssembleConcentrationRHS( UserContext *uc )
{
  int i;
  PetscReal *b, zeros[5] = {0,0,0,0,1};
  PetscInt bcidx, bcnodeidx;
  unsigned char color;
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  PetscLogEventBegin(EVENT_AssembleConcentrationRHS,0,0,0,0);
  
  VecSet(uc->b, 0.);
  VecGetArray(uc->b, &b);

  for (i = 0; i < uc->numBC; ++i)
  {
    bcidx = uc->bcToImage[i];
    color = uc->filedata[bcidx];
    bcnodeidx = uc->imageToNode[bcidx];
    b[bcnodeidx] = uc->BCLabels[color].concentrationBC;
    MatSetValues(uc->A, 1,&bcnodeidx,5,uc->nodes[bcnodeidx].star,zeros,INSERT_VALUES);
  }
  
  VecRestoreArray(uc->b, &b);
  
  ierr = MatAssemblyBegin(uc->A, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(uc->A, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  
  PetscLogEventEnd(EVENT_AssembleConcentrationRHS,0,0,0,0);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "SolveConcentration"
PetscInt EVENT_SolveConcentration;
PetscErrorCode SolveConcentration( UserContext *uc )
{
  PetscInt maxIndex;
  PetscReal vx, vy, CFL, dt, vel;
  PetscErrorCode ierr;
    
  PetscFunctionBegin;
  PetscLogEventBegin(EVENT_SolveConcentration,0,0,0,0);
  ierr = VecMax(uc->ss, &maxIndex, &vx); CHKERRQ(ierr);
  ierr = VecGetValues(uc->u, 1, &maxIndex, &vx); CHKERRQ(ierr);
  ierr = VecGetValues(uc->v, 1, &maxIndex, &vy); CHKERRQ(ierr);
  
  vel = PetscAbs(vx) + PetscAbs(vy);
  
  CFL = vel * DELTA_T / DELTA_X;
  
  dt  = DELTA_X / vel;  
  
  printf("CURRENT CFL: \t%f\n", CFL);
  printf("CURRENT dt:  \t%f\n", DELTA_T);
  printf("SUGGESTED dt:\t%f\n", dt);
  
  PetscLogEventEnd(EVENT_SolveConcentration,0,0,0,0);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "AssembleSSConcentration"
PetscInt EVENT_AssembleSSConcentration;
PetscErrorCode AssembleSSConcentration(UserContext *uc)
{
  int i, j;
  PetscReal *u, *v, ut, vt, val[5];
  const PetscReal Dh2 = DIFFUSION  / ( DELTA_X * DELTA_X );
  Node *n;
  PetscErrorCode ierr;;
  
  PetscFunctionBegin;
  PetscLogEventBegin(EVENT_AssembleSSConcentration,0,0,0,0);
  PetscLogEventRegister(&EVENT_AssembleSSConcentration,"AssembleSSConcentration", 0);
  VecGetArray(uc->u, &u);
  VecGetArray(uc->v, &v);
  
  for (i = 0; i < uc->numNodes; ++i)
  {
    n = &uc->nodes[i];

    // Diffusive part
    for( j = 0; j < 4; j++)
      val[j] = -1 * Dh2;
    val[4] = n->numNei * Dh2;

    // Advective part (upwinding) 
    ut = u[i] / DELTA_X;
    vt = v[i] / DELTA_X;
    
    if( u[i] < 0. ) { val[2] += ut; val[4] -= ut; }
    if( u[i] > 0. ) { val[1] -= ut; val[4] += ut; }
    if( v[i] < 0. ) { val[3] += vt; val[4] -= vt; }
    if( v[i] > 0. ) { val[0] -= vt; val[4] += vt; }

    ierr = MatSetValues(uc->A, 1, &i, 5, n->star, val, INSERT_VALUES); CHKERRQ(ierr);
  }
  
  VecRestoreArray(uc->v, &v);
  VecRestoreArray(uc->u, &u);
  PetscLogEventEnd(EVENT_AssembleSSConcentration,0,0,0,0);
  PetscFunctionReturn(0);
}
