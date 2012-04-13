#include "DWorld.h"
#include "petscts.h"
#include "MyCheck.h"



PetscErrorCode InitialCondition( DA da, Vec vec_chem );

START_TEST( test_DWorldTimeStepping )
{
  DWorld world;
  TS ts;
  PetscInt max_steps = 2;
  PetscReal max_time = 1500.0;
  PetscReal dt   = 1;
  PetscErrorCode ierr;
  
  Reaction rxnEC;
  ReactionCreate_GrayScott( &rxnEC );
    
  DWorldCreate(rxnEC, 25,25,3, 2.5,2.5,0.01, &world);
  InitialCondition_GrayScott(world->da, world->global);  

  
  ierr = TSCreate(PETSC_COMM_WORLD,&ts);CHKERRQ(ierr);
  ierr = TSSetProblemType(ts,TS_NONLINEAR);CHKERRQ(ierr);
  ierr = TSSetType(ts,TS_BEULER);CHKERRQ(ierr);
  ierr = TSSetInitialTimeStep(ts,0.0,dt);CHKERRQ(ierr);
  ierr = TSSetDuration(ts,max_steps,max_time);CHKERRQ(ierr);
  ierr = TSMonitorSet(ts,Monitor,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr); 
  ierr = TSSetFromOptions(ts);CHKERRQ(ierr);
  ierr = TSSetRHSFunction(ts,RHSFunction,world);CHKERRQ(ierr);
  ierr = TSSetRHSJacobian(ts,world->J,world->J,RHSJacobian,world);CHKERRQ(ierr);
  ierr = TSSetSolution(ts,world->global);CHKERRQ(ierr);
  
  SNES snes;
  ierr = TSGetSNES(ts, &snes); CHKERRQ(ierr);
  SNESSetFromOptions(snes);
  KSP ksp; 
  ierr = SNESGetKSP(snes,&ksp); CHKERRQ(ierr);
  KSPSetFromOptions(ksp);
  PC pc;
  ierr = KSPGetPC(ksp, &pc); CHKERRQ(ierr);
  
  KSPSetType(ksp, KSPPREONLY);
  PCSetType(pc,PCLU);
  
//  KSPView(ksp, PETSC_VIEWER_STDOUT_WORLD);
  
  ierr = TSSetUp(ts);CHKERRQ(ierr);
  
 
/* 
  KSPSetUp(ksp);
  
mark_point();
  KSP *subksp;
 
  PCBJacobiGetSubKSP(pc,0,0,&subksp);
 
  KSPGetPC(subksp[0],&pc);
  PCSetType(pc,PCLU);
  KSPSetType(subksp[0],KSPPREONLY);
  Mat pmat;
  KSPGetOperators(subksp[0],0,&pmat, 0);
  Mat dmat;
  MatConvert(pmat,MATUMFPACK,MAT_INITIAL_MATRIX,&dmat);
  KSPSetOperators(subksp[0],dmat,dmat,DIFFERENT_NONZERO_PATTERN);
*/
//  PC pc;
//  ierr = TSSundialsGetPC(ts,&pc);CHKERRQ(ierr);
//  ierr = PCSetType(pc,PCPBJACOBI);CHKERRQ(ierr);

//  ierr = TSSetUp(ts);CHKERRQ(ierr);
  
  PetscInt step = max_steps;
  PetscReal final_time = max_time;
  ierr = TSStep(ts,&step,&final_time);CHKERRQ(ierr);
  
  
}
END_TEST

START_TEST( view_Jacobian )
{
  int num=2;
  Reaction rxn;
  DWorld world;
  MatStructure flag;
  PetscErrorCode ierr;
  
  ReactionCreate_SimpleDegradation(&rxn);
  
  DWorldCreate(rxn,3,3,3,2,2,2,&world);
  
  RHSJacobian(0,0,world->global,&world->J,0,&flag,world);
  
  PetscViewer view;
  PetscViewerASCIIOpen(PETSC_COMM_SELF, 
      "/home/abergman/Research/DCell/temp/mat.dat", &view);
  MatView(world->J, view);  
  
  DWorldDestroy(world);
}
END_TEST

START_TEST( view_InitialCondition )
{
  PetscViewer binv;
  Reaction rxn;
  DWorld world;
  PetscErrorCode ierr;
 
  ReactionCreate(2, &rxn);
  DWorldCreate(rxn, 25, 25, 3, 2.5, 2.5, 1, &world);
  mark_point();
  InitialCondition_GrayScott(world->da, world->global);
  
  PetscViewerBinaryOpen(PETSC_COMM_SELF,
      "/home/abergman/Research/DCell/temp/init.dat",FILE_MODE_WRITE,&binv);
  VecView(world->global, binv);
  DWorldDestroy(world);
}
END_TEST

START_TEST( test_DAVecGetArrayDOF )
{
  int i, j, k, c;
  PetscInt xs,ys,zs,xm,ym,zm,mx,my,dof=3;
  PetscReal *array, ****grid;
  PetscInt len;
  Reaction rxn;
  DWorld world;
  PetscErrorCode ierr;
  int size;
  MPI_Comm_size(PETSC_COMM_WORLD, &size);
  if( size != 1 ) return;
  ReactionCreate( dof, &rxn);
  DWorldCreate(rxn, 4,5,1,1,1,0,&world);
  VecGetSize(world->global, &len);
  VecGetArray(world->global, &array );
  for (int i = 0; i < len; ++i) {
    array[i] = i;
  }
  VecRestoreArray(world->global, &array);
  DAGetCorners(world->da,&xs,&ys,&zs,&xm,&ym,&zm);
  DAVecGetArrayDOF(world->da, world->global, &grid);
  int count = 0;
  for( k=zs; k<zs+zm; ++k)
  {
    for (j=ys; j<ys+ym; j++) 
    {
      for (i=xs; i<xs+xm; i++)
      {
        for (c = 0; c < dof; ++c) {
          fail_unless( grid[k][j][i][c] == count, 
              "Grid indexing wrong: %f vs %d", grid[k][j][i][c], count);
          count++;
        }
      }
    }    
  }
  
  DAVecRestoreArrayDOF(world->da, world->global, &grid);
  DWorldDestroy(world);
}
END_TEST

START_TEST( test_DWorldCreateDestroy )
{
  Reaction rxn;
  DWorld world;
  PetscErrorCode ierr;
  
  ReactionCreate( 3, &rxn);
  DWorldCreate(rxn, 10,10,1,1,1,0,&world);
  DWorldDestroy(world);
}
END_TEST

Suite *DWorld_suite()
{
  Suite *s = suite_create ("DWorld Check");
  TCase *tc_core = tcase_create("Core");
  tcase_add_test( tc_core,  test_DWorldCreateDestroy );
  tcase_add_test( tc_core,  test_DAVecGetArrayDOF );
  tcase_add_test( tc_core,  view_InitialCondition );
  tcase_add_test( tc_core,  view_Jacobian );
  tcase_add_test( tc_core,  test_DWorldTimeStepping );
//  tcase_add_test( tc_core,  );
  tcase_set_timeout(tc_core, 0);
  suite_add_tcase( s, tc_core);
  return s;
}


PetscErrorCode InitialCondition( DA da, Vec vec_chem )
{
  PetscReal ****chem;
  PetscInt xs,ys,zs, xm,ym,zm,mx,my,dof;
  int i,j,k,c;
  PetscErrorCode ierr;
  
  ierr = DAGetInfo(da,0,&mx,&my,0,0,0,0,&dof,0,0,0);CHKERRQ(ierr);
  ierr = DAGetCorners(da,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
  ierr = DAVecGetArrayDOF(da,vec_chem,&chem);CHKERRQ(ierr);
  
  for (k=zs; k<zs+zm; ++k)
  {
    for (j=ys; j<ys+ym; j++) 
    {
      for (i=xs; i<xs+xm; i++)
      {
        for( c = 0; c < dof; ++c)
        {
          if( 0.25*mx < i && i < 0.75*mx && 0.25*my < j && j < 0.75*my )
            chem[k][j][i][c] = 1;
          else
            chem[k][j][i][c] = 0;
        }
      }
    }
  }
  
  ierr = DAVecRestoreArrayDOF(da,vec_chem,&chem);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}