#include "Main.h"
#include "MyCheck.h"
#include "PoissonSolver.h"

//#include "ImmersedInterfaceMethod.h"

START_TEST( test_Modulus )
{
  int d1=7;
  PetscErrorCode ierr;
  
  for (int i = 0; i < d1; ++i) 
  {
    printf("%d: %d, %d\n", i, (i+1)%d1, (i-1+d1)%d1);
  }
}
END_TEST

START_TEST( test_PeriodicBC )
{
  PetscErrorCode ierr;
  PetscInt d1 = 7, d2 = d1;
  Mat m;
  PetscViewer view;
  
  Generate2DLaplacianPeriodicBC( d1, d2, &m);
  
  ierr = MatView( m, PETSC_VIEWER_STDOUT_SELF); CHKERRQ(ierr);
  PetscViewerASCIIOpen(PETSC_COMM_SELF, "mat.dat", &view);
  MatView(m, view);
  
  MatDestroy(m);
}
END_TEST
/*
START_TEST( test_SpectralIIM )
{
  PetscErrorCode ierr;
  int d1 = 211, d2 = d1;
  LevelSet2D ls;
  CreateLevelSet2D(d1, d2, &ls);
  LevelSetInitializeToStar(ls,PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE);
  WriteVector("initls",ls->g2d->v);
  Grid2D rhsC;
  CreateGrid2D(d1, d2, &rhsC);
  Grid2D px, py;
  CreateGrid2D(d1, d2, &px);
  CreateGrid2D(d1, d2, &py);
LevelSet2D lsnew;
CreateLevelSet2D(d1, d2, &lsnew);
WriteVec wv, wvx,wvy;
WriteVecCreate(lsnew->g2d->v, "ls",&wv);
WriteVecCreate(px->v, "vx",&wvx);
WriteVecCreate(py->v, "vy",&wvy);

for( int t=0; t<100; t++) {
  printf("time: %d\n",t);
  ComputeCorrection(ls, rhsC->v);
  
  WriteVector("rhscorrec",rhsC->v);
  SpectralMethod2D( rhsC );
  WriteVector("sol",rhsC->v);
  
  int i, j;
  for( j = 1; j < d2-1; ++j)
  {
    for( i = 1; i < d1-1; ++i)
    {
      px->v2[i][j] = (rhsC->v2[i+1][j] - rhsC->v2[i-1][j]) / 2.;
      py->v2[i][j] = (rhsC->v2[i][j+1] - rhsC->v2[i][j-1]) / 2.;
    }
  }
  mark_point();
  PetscReal **phi = ls->g2d->v2;
    IrregularNode n;
  for( i=0; i<ls->irregularNodes->len; i++)
  {
    n = g_array_index( ls->irregularNodes, IrregularNode, i);
    if( n.sign == PetscSign( phi[n.x+1][n.y] ) ) {
      px->v2[n.x][n.y] = rhsC->v2[n.x+1][n.y] - rhsC->v2[n.x][n.y];
    } else if( n.sign == PetscSign( phi[n.x-1][n.y] ) ) {
      px->v2[n.x][n.y] = rhsC->v2[n.x][n.y] - rhsC->v2[n.x-1][n.y];
    } else {
//      printf("Opposite sign on both sides of node\n");
//      exit(1);
    }
    if( n.sign == PetscSign( phi[n.x][n.y+1] ) ) {
      py->v2[n.x][n.y] = rhsC->v2[n.x][n.y+1] - rhsC->v2[n.x][n.y];
    } else if( n.sign == PetscSign( phi[n.x][n.y-1] ) ) {
      py->v2[n.x][n.y] = rhsC->v2[n.x][n.y] - rhsC->v2[n.x][n.y-1];
    } else {
//      printf("Opposite sign on both sides of node\n");
//      exit(1);
    }
  }
  mark_point();
//  WriteVector("px",px->v);
//  WriteVector("py",py->v);
  PetscReal imu=1000;
  VecScale(px->v, imu);
  VecScale(py->v, imu);
  SpectralMethod2D( px );
  SpectralMethod2D( py );
  WriteVecToDisk(wvx);
  WriteVecToDisk(wvy);
  
  PetscReal maxX,maxY;
  VecNorm(px->v, NORM_INFINITY, &maxX);
  VecNorm(py->v, NORM_INFINITY, &maxY);
  LevelSetAdvect2D(1/(maxX+maxY), px, py, ls, lsnew);
  WriteVecToDisk(wv);
  VecCopy( lsnew->g2d->v, ls->g2d->v);
  UpdateIrregularNodeList(ls);
  ReinitializeLevelSet(ls);
}
}
END_TEST
*/
START_TEST( viz_Dirichlet )
{
   PetscErrorCode ierr;
   SpectralMethod2DTimer();
}
END_TEST
START_TEST( test_SpectralMethod )
{
  PetscErrorCode ierr;
  PetscLogDouble t;
  for( int i = 4; i < 256; i+=2 )
  {
    ierr = SpectralMethod(i, &t); CHKERRQ(ierr);
    printf("{ %d, %f},\n",i,t);
  }
}
END_TEST

START_TEST( GenerateLaplacian2DNoBC_test )
{
  PetscInt d1, d2;
  Mat m;
  PetscViewer view;
  
  d1 = d2 = 1e2;
  
  GenerateLaplacian2DNoBC(d1, d2, &m);
  PetscViewerASCIIOpen(PETSC_COMM_SELF, "mat.dat", &view);
  MatView(m, view);
}
END_TEST

START_TEST( Generate2DLapacian_test )
{
  PetscErrorCode ierr;
  PetscInt d = 50, len = d*d;
  Mat m;
  PetscViewer view;
  Generate2DLapacian( d, d, &m);
  
  PetscViewerASCIIOpen(PETSC_COMM_SELF, "mat.dat", &view);
  MatView(m, view);
  
  KSP ksp;
  PC pc;
  ierr = KSPCreate(PETSC_COMM_SELF, &ksp); CHKERRQ(ierr);
  KSPSetType(ksp, KSPPREONLY);
  KSPGetPC(ksp, &pc);
  PCSetType(pc, PCCHOLESKY);
  KSPSetOperators(ksp, m, m, SAME_PRECONDITIONER);
  Vec v, s;
  VecCreateSeq(PETSC_COMM_SELF,len,&v);
  VecDuplicate(v, &s);
  VecSet(v, 1);
  ierr = KSPSolve(ksp, v, s); CHKERRQ(ierr);
  
}
END_TEST

Suite* Poisson_suite (void)
{
  Suite *s = suite_create ("Poisson Solver Check");
  TCase *tc_core = tcase_create("Core");
//  tcase_add_test( tc_core,  Generate2DLapacian_test );
//  tcase_add_test( tc_core,  GenerateLaplacian2DNoBC_test );
//  tcase_add_test( tc_core,  test_SpectralMethod );
  tcase_add_test( tc_core,  viz_Dirichlet );
//  tcase_add_test( tc_core,  test_SpectralIIM );
//  tcase_add_test( tc_core,  test_PeriodicBC );
//  tcase_add_test( tc_core,  test_Modulus );
  tcase_set_timeout(tc_core, 0);
  suite_add_tcase( s, tc_core);
  return s;
}

Suite* IntelFFT_Suite(void);

int RunCheck()
{ 
  Suite *s = Poisson_suite ();
  SRunner *sr = srunner_create (IntelFFT_Suite());
//  srunner_add_suite(sr, IntelFFT_Suite() );
  srunner_set_fork_status(sr, CK_NOFORK);
  srunner_run_all (sr, CK_NORMAL);
  srunner_free (sr);
  return 0;
}

PetscErrorCode PetscMain()
{
  PetscFunctionReturn(0);
}


START_TEST( test_Diffusion )
{
  int nn = 128;
  IntelFFT fft;
  double D = 1, dt = .001, q = 1 / (D*dt), X = 1;
  Vec F, source, copy;
  PetscErrorCode ierr;
  
  ierr = IntelFFTCreate(nn,nn,nn,&fft, &F); CHKERRQ(ierr);
  ierr = IntelFFTSetDomain(fft,0,X,0,X,0,X); CHKERRQ(ierr);
  IntelFFTSetType(fft, q, 'N');
  VecDuplicate(F,&source);
  VecDuplicate(F,&copy);
  int half = (nn)/2;
  int c = half + half * (nn+1) + half * (nn+1)*(nn+1);
  printf("c:%d\n",c);
  double v = nn*nn;
  VecSetValues(source,1,&c,&v,INSERT_VALUES);
//  VecSet(source,1.);
  
  /*
  for( int i = 0; i < nn / 4; i++ )
    for( int j = 0; j < nn / 4; j++ )
      for( int k = 0; k < nn / 4; ++k)
      {
        c = i + j * (nn+1) + k * (nn+1)*(nn+1);
        VecSetValues(F,1,&c,&v,INSERT_VALUES);
      }
  */
  char tmp[128], format[128];
  PetscGetTmp(PETSC_COMM_SELF,tmp,128 );
  
  PetscLogDouble p1, p2, pToT = 0;
  PetscGetTime(&p1);
  int steps = 10;
  for( int i = 0; i < steps; i++ )
  {
    VecAYPX(F,q,source);
    
//    VecScale(F,q);
    
    
    ierr = IntelFFTSolve(fft); CHKERRQ(ierr);
    
    
//    VecCopy(F,copy);
//    VecCopy(F,copy);
    
    
    sprintf(format,"%s/test.%d", tmp,i);
    
    PetscViewer binv;
    PetscViewerBinaryOpen(PETSC_COMM_WORLD,format,FILE_MODE_WRITE,&binv);
    VecView(F, binv);
    PetscViewerDestroy(binv);
  }
  PetscGetTime(&p2);
  printf("TIME: %f\n", (p2-p1));
  pToT += (p2-p1);
  
  printf("TOT: %f\tAVG: %f\n", pToT, pToT / steps);
  
  
  
  ierr = IntelFFTDestroy(fft); CHKERRQ(ierr);
}
END_TEST

START_TEST( test_CreateDestroy )
{
  int nn = 32;
  Vec F;
  IntelFFT fft;
  PetscErrorCode ierr;

  ierr = IntelFFTCreate(nn,nn,nn,&fft, &F); CHKERRQ(ierr);
  ierr = IntelFFTDestroy(fft); CHKERRQ(ierr);
}
END_TEST

Suite* IntelFFT_Suite(void)
{
  Suite *s = suite_create ("IntelFFT Suite");
  TCase *tc_core = tcase_create("Core");
  tcase_add_test( tc_core, test_CreateDestroy );
  tcase_add_test( tc_core, test_Diffusion );
  //  tcase_add_test( tc_core,  );  
  tcase_set_timeout(tc_core, 0);
  suite_add_tcase( s, tc_core);

  return s;
}
