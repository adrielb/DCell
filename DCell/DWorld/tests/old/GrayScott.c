#include "DCell.h"


int main(int argc, char **args)
{
  PetscErrorCode  ierr;
  ierr = PetscInitialize(&argc, &args, (char *) 0, ""); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "Start\n"); CHKERRQ(ierr);
  
  DWorld world;
  TS ts;
  PetscInt step, max_steps = 15;
  PetscReal time, max_time = 15000.0;
  PetscReal dt   = 1;
  
  Reaction rxnEC;
  ReactionCreate_GrayScott( &rxnEC );
  /* 3D
  iCoor num = {100,100,25};
  Coor len  = {2.5, 2.5, 0.625};
  */
  iCoor num = {100,100,0};
  Coor len  = {2.5, 2.5,0};
  
  DWorldCreate(rxnEC, num, len, &world);
  InitialCondition_GrayScott2D(world->da, world->global);
  
  /* Load check-point
  PetscViewer in;
  PetscViewerBinaryOpen(PETSC_COMM_WORLD,"/home/abergman/Research/DCell/temp/vec.2220.dat",FILE_MODE_READ,&in);
  VecLoadIntoVector(in,world->global);
  PetscViewerDestroy(in);
  */
  
  ierr = TSCreate(PETSC_COMM_WORLD,&ts);CHKERRQ(ierr);
  ierr = TSSetProblemType(ts,TS_NONLINEAR);CHKERRQ(ierr);
//  ierr = TSSetType(ts,TS_CRANK_NICHOLSON);CHKERRQ(ierr);
  ierr = TSSetType(ts,TS_BEULER);CHKERRQ(ierr);
  ierr = TSSetRHSFunction(ts,RHSFunction,world);CHKERRQ(ierr);
  ierr = TSSetRHSJacobian(ts,world->J,world->J,RHSJacobian,world);CHKERRQ(ierr);
  ierr = TSSetInitialTimeStep(ts,0.0,dt);CHKERRQ(ierr);
  ierr = TSSetDuration(ts,max_steps,max_time);CHKERRQ(ierr);
  ierr = TSMonitorSet(ts,Monitor,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr); 
  ierr = TSSetFromOptions(ts);CHKERRQ(ierr);
  
  ierr = TSSetSolution(ts,world->global);CHKERRQ(ierr);
  ierr = TSSetUp(ts);CHKERRQ(ierr);
  
  SNES snes;
  ierr = TSGetSNES(ts, &snes); CHKERRQ(ierr);
  SNESSetFromOptions(snes);
  KSP ksp; 
  ierr = SNESGetKSP(snes,&ksp); CHKERRQ(ierr);
  KSPSetFromOptions(ksp);
  PC pc;
  ierr = KSPGetPC(ksp, &pc); CHKERRQ(ierr);
  
//  KSPSetType(ksp, KSPPREONLY);
//  PCSetType(pc,PCLU);
  
  KSPSetType(ksp, KSPGMRES);
  PCSetType(pc, PCBJACOBI);
  
  ierr = TSSetDuration(ts,1,1);CHKERRQ(ierr);
  ierr = TSStep(ts,&step,&time);CHKERRQ(ierr);
  ierr = TSSetDuration(ts,max_steps,max_time);CHKERRQ(ierr);
  
/*
 * Using an iterative solver on each block
 */
  KSP *subksp;
  PCBJacobiGetSubKSP(pc,0,0,&subksp);
  KSPSetType(subksp[0],KSPGMRES);
  KSPGetPC(subksp[0],&pc);
  PCSetType(pc,PCILU); //TODO: try a pc type of SOR
  
  KSPView(ksp, PETSC_VIEWER_STDOUT_WORLD);
  printf("====================");
  KSPView(subksp[0], PETSC_VIEWER_STDOUT_WORLD);
  

/* using a Direct Solver 
 * 
 * KSPSetUp(ksp);
    
  KSP *subksp;

  PCBJacobiGetSubKSP(pc,0,0,&subksp);
  KSPView(subksp[0], PETSC_VIEWER_STDOUT_WORLD);
  KSPSetType(subksp[0],KSPPREONLY);
  KSPGetPC(subksp[0],&pc);
  PCSetType(pc,PCLU);
  Mat pmat;
  KSPGetOperators(subksp[0],0,&pmat, 0);
  Mat dmat;
  MatConvert(pmat,MATUMFPACK,MAT_INITIAL_MATRIX,&dmat);
  KSPSetOperators(subksp[0],dmat,dmat,DIFFERENT_NONZERO_PATTERN);
  
*/
  step = max_steps;
  time = max_time;
  ierr = TSStep(ts,&step,&time);CHKERRQ(ierr);
  
  ierr = DWorldDestroy(world); CHKERRQ(ierr);

  ierr = PetscPrintf(PETSC_COMM_WORLD, "End\n"); CHKERRQ(ierr);
  ierr = PetscFinalize(); CHKERRQ(ierr);
  return 0;
}
PetscErrorCode InitialCondition_GrayScott2D( DA da, Vec vec_chem )
{
  PetscReal X, Y;
  PetscReal ***chem;
  PetscInt xs,ys,zs, xm,ym,zm,mx,my,mz,dof;
  int i,j,c;
  PetscErrorCode ierr;
  
  ierr = DAGetInfo(da,0,&mx,&my,&mz,0,0,0,&dof,0,0,0);CHKERRQ(ierr);
  ierr = DAGetCorners(da,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
  ierr = DAVecGetArrayDOF(da,vec_chem,&chem);CHKERRQ(ierr);
  
  
  for (j=ys; j<ys+ym; j++) 
  {
    for (i=xs; i<xs+xm; i++)
    {
      X = i * 2.5 / (mx-2);
      Y = j * 2.5 / (my-2);
      
      // V = 0.25 * Sin( 4 * pi * X)^2 * Sin(4 * pi * Y)^2 
      if( 1 < X && X < 1.5 && 1 < Y && Y < 1.5) // 1 < X,Y < 1.5
        chem[j][i][1] = 0.25 * 
          PetscSqr(sin(4*PETSC_PI*X)) * 
          PetscSqr(sin(4*PETSC_PI*Y));
      else 
        chem[j][i][1] = 0;
      
      // U = 1 - 2 * V
      chem[j][i][0] = 1 - 2 * chem[j][i][1]; 
    }
  }
  
 
  ierr = DAVecRestoreArrayDOF(da,vec_chem,&chem);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
PetscErrorCode InitialCondition_GrayScott( DA da, Vec vec_chem )
{
  PetscReal X, Y;
  PetscReal ****chem;
  PetscInt xs,ys,zs, xm,ym,zm,mx,my,mz,dof;
  int i,j,k,c;
  PetscErrorCode ierr;
  
  ierr = DAGetInfo(da,0,&mx,&my,&mz,0,0,0,&dof,0,0,0);CHKERRQ(ierr);
  ierr = DAGetCorners(da,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
  ierr = DAVecGetArrayDOF(da,vec_chem,&chem);CHKERRQ(ierr);
  
  for( k=zs; k<zs+zm; ++k)
  {
    for (j=ys; j<ys+ym; j++) 
    {
      for (i=xs; i<xs+xm; i++)
      {
        X = i * 2.5 / (mx-1);
        Y = j * 2.5 / (my-1);
        
        // V = 0.25 * Sin( 4 * pi * X)^2 * Sin(4 * pi * Y)^2 
        if( 1 < X && X < 1.5 && 1 < Y && Y < 1.5) // 1 < X,Y < 1.5
          chem[k][j][i][1] = 0.25 * 
            PetscSqr(sin(4*PETSC_PI*X)) * 
            PetscSqr(sin(4*PETSC_PI*Y)) * 
            -k * (k - (mz-1)) / (PetscSqr(mz-1)/4);
        else 
          chem[k][j][i][1] = 0;
        
        // U = 1 - 2 * V
        chem[k][j][i][0] = 1 - 2 * chem[k][j][i][1]; 
      }
    }
  }
 
  ierr = DAVecRestoreArrayDOF(da,vec_chem,&chem);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
