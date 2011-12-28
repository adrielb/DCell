#include "MyCheck.h"
#include "DCell.h"

START_TEST( viz_AssembleDiffusion )
{
  int d1=64, d2=d1, d3=d2;
  PetscErrorCode ierr;
  LevelSet3D ls;
  LevelSet3DCreate("ls",d1,d2,d3,&ls);
  LevelSetInitializeToSphere(ls, PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE);
  WriteVecToDisk(ls->g3d->wv);
  
  DA da;
  DACreate3d(PETSC_COMM_SELF,//MPI Communicator   
    DA_NONPERIODIC,   //DA_NONPERIODIC, DA_XPERIODIC, DA_YPERIODIC, DA_XYPERIODIC
    DA_STENCIL_STAR, //DA_STENCIL_BOX or DA_STENCIL_STAR
    d1, d2, d3,     //Global dimension
    1, 1, 1,       //Number procs per dim
    1,            //dof
    1,           //stencil width
    0,0,0,      //specific array of nodes
    &da);
  
  Mat mat;
  DAGetMatrix( da, MATSEQAIJ, &mat);
  
  DCellAssembleDiffusion(ls, mat);

/*
  PetscViewer view;
  PetscViewerASCIIOpen(PETSC_COMM_SELF, "/home/abergman/Research/DCell/temp/mat.dat", &view);
  MatView(mat, view);
  
  Grid3D diag;
  Grid3DCreate("diag",d1,d2,d3,&diag);
  MatGetDiagonal(mat, diag->v);
  WriteVecToDisk(diag->wv);
*/

  PetscTruth truth;
  MatIsSymmetric(mat, .001, &truth);
  fail_unless(truth, "mat not symmetric!");
  
  KSP ksp;
  KSPCreate(PETSC_COMM_SELF,&ksp);
  PC pc;
  KSPGetPC(ksp, &pc);
  KSPSetOperators(ksp,mat, mat, SAME_PRECONDITIONER );
  PCFactorSetMatOrderingType(pc, MATORDERING_ND);
  
  KSPSetType(ksp, KSPCG);
  PCSetType(pc, PCICC);
  PCFactorSetLevels(pc, 4);
  KSPSetInitialGuessNonzero(ksp, PETSC_TRUE);
  
//  KSPSetType(ksp,KSPPREONLY);
//  PCSetType(pc, PCCHOLESKY);

  Grid3D c;
  Grid3DCreate("chem", d1,d2,d3, &c );
  VecZeroEntries(c->v);
  c->v3[d1/2][d2/2][d3/2] = d1*d2*d3;
  WriteVecToDisk(c->wv);    

  for (int i = 0; i < 333; ++i)
  {
    for( int j = 0; j < 1; ++j)
    {
      KSPSolve(ksp, c->v, c->v);      
    }
    WriteVecToDisk(c->wv);    
  }
  
  
  IrregularNodeListWrite(ls->irregularNodes,0);
}
END_TEST

Suite *Reaction_suite();
Suite *DWorld_suite();
Suite* CreateSuite (void)
{
  Suite *s = suite_create ("DCell Check");

  TCase *tc_core = tcase_create("Core");
//  tcase_add_test( tc_core,  viz_AssembleDiffusion );
 // tcase_add_test( tc_core,   );
  suite_add_tcase( s, tc_core);  
  return s;
}

int RunCheck()
{ 
  Suite *s = CreateSuite ();
  SRunner *sr = srunner_create (s);
  srunner_add_suite(sr, Reaction_suite() );
  srunner_add_suite(sr, DWorld_suite() );
  srunner_set_fork_status(sr, CK_NOFORK);
  srunner_run_all (sr, CK_NORMAL);
  srunner_free (sr);
  return 0;
}