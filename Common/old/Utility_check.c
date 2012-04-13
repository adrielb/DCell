#include "mkl.h"
#include "Main.h"
#include "MyCheck.h"
#include "Utilities.h"
#include "Grid.h"
#include "petscksp.h"
#include "petscda.h"

double identity3d(double ***p, int x, int y, int z )
{
  return p[x][y][z];
}

START_TEST( test_Bilinear3D )
{
  PetscErrorCode ierr;
  int i,j,k,c=0;
  int d1=3,d2=d1,d3=d1;
  int f = 3;
  Grid3D g,gb;
  ierr = Grid3DCreate("test",d1,d2,d3,&g); CHKERRQ(ierr);
  ierr = Grid3DCreate("bi",d1*f,d2*f,d3*f,&gb); CHKERRQ(ierr);
  for (i = 0; i < g->len; ++i)
  {
    g->v1[i] = i;
  }
  PetscReal dh = (d1 - 1.01 - .01)/(f*d1 - 1 - 0);
  PetscReal  a = 0.001;
  for (i = 0; i < gb->d1; ++i)
  {
    for (j = 0; j < gb->d2; ++j)
    {
      for (k = 0; k < gb->d3; ++k)
      {
        gb->v3[i][j][k] = Bilinear3D(identity3d, g, i*dh+a,j*dh+a,k*dh+a);
      }
    }
  }
  WriteVecToDisk(gb->wv);
}
END_TEST

START_TEST( test_Grid3D )
{
  PetscErrorCode ierr;
  int i,j,k,c=0;
  int d1=3,d2=4,d3=5;
  Grid3D g;
  ierr = Grid3DCreate("test",d1,d2,d3,&g); CHKERRQ(ierr);
  for (i = 0; i < g->len; ++i)
  {
    g->v1[i] = i;
  }

  DA da;
  ierr = DACreate3d(PETSC_COMM_SELF,//MPI Communicator   
    DA_NONPERIODIC,   //DA_NONPERIODIC, DA_XPERIODIC, DA_YPERIODIC, DA_XYPERIODIC
    DA_STENCIL_STAR, //DA_STENCIL_BOX or DA_STENCIL_STAR
    d1, d2, d3,     //Global dimension
    1, 1, 1,       //Number procs per dim
    1,            //dof
    1,           //stencil width
    0,0,0,      //specific array of nodes
    &da);
  
  Mat m;  
  ierr = DAGetMatrix( da, MATSEQAIJ, &m); CHKERRQ(ierr);
  
  PetscReal val=0;
  MatStencil row;row.c=0;
  for (row.k = 0; row.k < g->d3; ++row.k)
  {
    for (row.j = 0; row.j < g->d2; ++row.j)
    {
      for (row.i = 0; row.i < g->d1; ++row.i)
      {
        ierr = MatSetValuesStencil(m, 1, &row, 1, &row, &val,INSERT_VALUES); CHKERRQ(ierr);
        val++;
      }
    }
  }
  MatAssemblyBegin(m,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(m,MAT_FINAL_ASSEMBLY);
//  MatView(m, PETSC_VIEWER_STDOUT_SELF);
  Grid3D diag;
  Grid3DCreate("diag",d1,d2,d3,&diag);
  MatGetDiagonal(m, diag->v);
//  VecView(diag->v, PETSC_VIEWER_STDOUT_SELF);
  
  c=0;
  for (i = 0; i < g->d1; ++i)
  {
    for (j = 0; j < g->d2; ++j)
    {
      for (k = 0; k < g->d3; ++k)
      {
        fail_unless(g->v3[i][j][k] == g->v1[c], "triple index failed: %f != %f at (%d, %d, %d) or %d\n", g->v3[i][j][k], g->v1[c], i,j,k,c);
        val = g->v1[From3Dto1D(g,i,j,k)];
        fail_unless(val == g->v1[c], "single index failed: %f != %f at (%d, %d, %d) or %d\n", val, g->v1[c], i,j,k,c);
        fail_unless(val == diag->v1[c], "single index failed: %f != %f at (%d, %d, %d) or %d\n", val, diag->v1[c], i,j,k,c); 
        c++;
      }
    }
  }
  WriteVecToDisk(g->wv);
  ierr = Grid3DDestroy(g); CHKERRQ(ierr);
}
END_TEST

double identity(double **p, int x, int y )
{
  return p[x][y];
}

START_TEST( test_GridFunction )
{
  PetscReal f=4;
  int i, j, d1=4,d2=4;
  Grid2D g2d, b;
  CreateGrid2D(d1,d2, &g2d);
  CreateGrid2D(d1*f-1, d2*f-1, &b);
  
  
  for( i = 0; i < d1*d2; i++ )
    g2d->v1[i] = i;
  
  for (j = 0; j < d2 * f - 1; ++j)
  {
    for (i = 0; i < d1 * f - 1; ++i)
    {
      b->v2[i][j] = Bilinear2D(identity, g2d, i/f ,j/f);
    }
  }
  
  WriteVector("bilineartest",b->v);
}
END_TEST

START_TEST( test_WriteVec )
{
  WriteVec wv;
  Vec v;
  PetscErrorCode ierr;
  VecCreateSeq( PETSC_COMM_SELF, 5, &v);
  WriteVecCreate( v, "testWriteVec", &wv);
  for( int i = 0; i < 5; i++)
  {
    printf("\n%d\n", i);
    VecSet( v, i);
    WriteVecToDisk(wv);
  }
  WriteVecDestroy( wv );
  VecDestroy(v);
}
END_TEST

START_TEST( test_PetscTemp )
{
  PetscErrorCode ierr;
  size_t len = 100;
  char* dir;
  PetscMalloc( len, &dir);
  ierr = PetscGetTmp( PETSC_COMM_WORLD, dir, len); CHKERRQ(ierr);
  printf("\n\nDIR: %s\n\n", dir);
  ierr = PetscFree(dir); CHKERRQ(ierr);
  
  PetscReal data[4] = {6,2,3,4};
  ierr = WriteVectorArray("tempdir",4,data); CHKERRQ(ierr);
  ierr = WriteVectorArray("tempdir1",4,data); CHKERRQ(ierr);
  ierr = WriteVectorArray("tempdir2",4,data); CHKERRQ(ierr);
  ierr = WriteVectorArray("tempdir3",4,data); CHKERRQ(ierr);
}
END_TEST

START_TEST( test_LeastSq )
{
  PetscErrorCode ierr;
  PetscReal x[7] = {-2.50725, -1.7072, -0.726423, 0.452286, 1.81034, 2.53149, 3.76245},
            y[7] = {1, 2, 5, 9, 7, 6, 5},
            sol[3] = {6.715669164916147,1.2967934985903051,-0.9888109988616558};
  int np = 7;
  PetscReal *s, *g;
  
  LeastSq ls;
  ierr = LeastSqCreate( 17, &ls); CHKERRQ(ierr);
  ierr = LeastSqSetNumPoints(ls, 7); CHKERRQ(ierr);
  ierr = LeastSqGetVecs(ls, &s, &g, PETSC_NULL); CHKERRQ(ierr);
  for (int i = 0; i < np; ++i) {
    s[i] = x[i];
    g[i] = y[i];
  }
  ierr = LeastSqSolve( ls ); CHKERRQ(ierr);
  for (int i = 0; i < 3; ++i) {
    fail_unless(PetscAbs(g[i]-sol[i])<.0001, "g: %f\t s:%f\t\n", g[i], sol[i] );
  }
  
  PetscReal xx[4] = {1, 4, 8, 10},
            yy[4] = {1, 5, 3, 6},
            ssol[3]= {0.7202797202797178, 0.8846153846153864, -0.09090909090909123};
  ierr = LeastSqSetNumPoints(ls, 4); CHKERRQ(ierr);
  for (int i = 0; i < np; ++i) {
    s[i] = xx[i];
    g[i] = yy[i];
  }
  ierr = LeastSqSolve( ls ); CHKERRQ(ierr);
  for (int i = 0; i < 3; ++i)
    fail_unless(PetscAbs(g[i]-ssol[i])<.0001, "%d\t g: %f\t s:%f\t\n", i, g[i], ssol[i] );
  
  ierr = LeastSqDestroy(ls); CHKERRQ(ierr);  
}
END_TEST

START_TEST( test_To1D )
{
  Grid2D g;
  PetscReal d1 = 5, d2 = 7;
  int count = 0;
  PetscErrorCode ierr;
  
  CreateGrid2D(5,7,&g);
  for( int i = 0; i < g->len; i++ )
    g->v1[i] = i;
  for( int j = 0; j < g->d2; j++ )
  {
    for( int i = 0; i < g->d1; i++ )
    {
      fail_unless( g->v1[To1D(g, i, j)] == count, 
                   "g: %f, c: %d", g->v1[To1D(g, i, j)], count )
      count++;
    }
  }
  
}
END_TEST

START_TEST( using_MKL_DROT )
{
  int i;
  PetscErrorCode ierr;
  
  double *x, *y, c = 0.9, s = 0.1;
  int incx = 1, incy = 1, n = 5;
  
  PetscMalloc( n*sizeof(double), &x);
  PetscMalloc( n*sizeof(double), &y);
  
  for( i = 0; i < n; ++i)
  {
    x[i] = i;
    y[i] = one;
  }
  
  for( i = 0; i < n; ++i)
  {
    printf("{%f, %f},\n",  x[i],  y[i]);
  }
  
  drot(&n, x, &incx, y, &incy, &c, &s);
  
  for( i = 0; i < n; ++i)
  {
    printf("{%f, %f},\n",  x[i],  y[i]);
  }
  
  PetscFree( x );
  PetscFree( y );
}
END_TEST

START_TEST( using_MKL_DGEMV )
{
//  y:= alpha*a*x + beta*y
  PetscErrorCode ierr;
  char trans = 'n';
  int m = 3, n = 4;
  
  double alpha = one;
  double *a;
  PetscMalloc( m*n*sizeof(double), &a);
  for( int i = 0; i < m*n; i++)
    a[i] = i+1;
  int lda = m;
  
  double *x;
  PetscMalloc( n*sizeof(double), &x);
//  for(int i = 0; i < n; i++)
    x[0] = one;
    x[1] = zero;
    x[2] = zero;
    x[3] = zero;
  int incx = 1;
  
  double beta = zero;
  double *y;
  PetscMalloc( n*sizeof(double), &y);
  for(int i = 0; i < n; i++)
    y[i] = zero;
  int incy = 1;
  
  dgemv(&trans,&m,&n,&alpha,a,&lda,x,&incx,&beta,y,&incy);
  
  for( int i = 0; i < n; i++ )
    printf("y: %f\n", y[i]);
  
  PetscFree(a);
  PetscFree(x);
  PetscFree(y);
}
END_TEST

START_TEST( using_MKL_DOT )
{
  PetscErrorCode ierr;
  
  int n = 5, incx = 1, incy = 1;
  double *x, *y, res;
  PetscMalloc(n*sizeof(double), &x);
  PetscMalloc(n*sizeof(double), &y);
  
  for( int i = 0; i < n; i++)
  {
    x[i] = i+3;
    y[i] = i+6;
  }
  
  res = ddot( &n, x, &incx, y, &incy);
  
  fail_unless( res == 210. , "res is %f", res);
  
  PetscFree(x);
  PetscFree(y);
}
END_TEST

START_TEST( using_MKL_SVD )
{
  PetscErrorCode ierr;

  int Np = 7;  
  int m = Np, n = 3;
  int lda = m;
  double *a;
  PetscMalloc(m*n*sizeof(double), &a);
  double aa[21] = {1., -2.50725, 6.28632, 1., -1.7072, 2.91453, 1., -0.726423, 0.527691, 1., \
0.452286, 0.204562, 1., 1.81034, 3.27735, 1., 2.53149, 6.40844, 1., 3.76245, \
14.1561};
  double at[21] = {1, 1, 1, 1, 1, 1, 1, -2.50725, -1.7072, -0.726423, 0.452286, 1.81034, \
      2.53149, 3.76245, 3.14316, 1.45727, 0.263846, 0.102281, 1.63868, 3.20422, \
      7.07805};
  double y[7]={1, 2, 5, 9, 7, 6, 5};
  for( int i = 0; i < m*n; i++ )
    a[i] = at[i];
  
  int ldb = Np;
  int nrhs = 1;
  double *b;
  PetscMalloc(Np*nrhs*sizeof(double), &b);
  PetscMemzero(b,Np*nrhs*sizeof(double));
//  double bb[9] = {one,  zero, zero,
//                  zero, one,  zero,
//                  zero, zero, one};
//  for( int i = 0; i < nrhs*m; i++)
//    b[i] = bb[i];

/* NEED B matrix to be Np identity!!! */
  for( int i = 0; i < Np*nrhs; i++)
    b[i] = y[i];
//    b[i+i*Np] = one;
    
  double rcond = -1;
  
  int lwork = 1000;
  double *work;
  PetscMalloc(lwork*sizeof(double), &work);
//  PetscMemzero(work, lwork*sizeof(double));
  
  int rank;
  double *s;
  PetscMalloc( m*sizeof(double), &s);
  int info;
  
  dgelss(
    &m /* num rows in A */,   //3 
    &n /* num cols in A */,   //7
    &nrhs /* num cols in B */,//3
    a /* matrix A for SVD */, 
    &lda /* first dim of A = max( 1, m )*/,  //3 
    b /* matrix B of RHS (become solution matrix x) */, 
    &ldb /* first dim of B = max( 1, m, n ) */, //7
    s /* matrix of s.v */, 
    &rcond /* relative condition s.v. treated as zero */, //-1 
    &rank /*  */, 
    work /*  */, 
    &lwork /* size of work */, 
    &info /*  */
  );
  
  printf("rank: %d\n", rank);
  printf("work[0]: %f\n", work[0] );
  printf("info: %d\n", info);
  
  for( int i = 0; i < n+0*Np*nrhs; i++)
    printf("b %d: %f\n", i, b[i]);
  printf("\n");
  for( int i = 0; i < 3; i++)
    printf("s %d: %f\n", i, s[i]);
  printf("\n");
  for( int i = 0; i < m*n; i++)
    printf("a %d: %f\n", i, a[i]);
  printf("\n");
}
END_TEST

START_TEST( check_FloatAsIntArrayIndex )
{
  int count = 0;
  double f, c;
  PetscErrorCode ierr;
  
  for(double i = 0.5; i < 10000; ++i) 
  {
    f = floor(i);
    c = ceil(i);
    fail_unless(count == (int)f, "floor(%f) failed %d", f, count);
    fail_unless(count+1==(int)c, "ceil(%f) failed %d", c, count);
    count++;
  }
}
END_TEST

START_TEST( CreateGrid2D_test )
{
  Grid2D g;
	int count = 0;
  
  PetscErrorCode ierr;
	 
  ierr = CreateGrid2D(8,8,&g); CHKERRQ(ierr);
  for( int i = 0; i < g->d1; ++i) 
  {
    for( int j = 0; j < g->d2; ++j) 
    {
      g->v2[i][j] = count;
      count++;
    }
  }
  
  for (int i = 0; i < g->len; ++i) 
  {
    fail_unless( PetscAbs(i-g->v1[i])<1e-9, "i:%d\t v:%f\n",i,g->v1[i]);
  }
  
  Bilinear2D(GridFunction2D_DerivX, g, 4.2, 4.3);
  
  ierr = DestroyGrid2D(g); CHKERRQ(ierr);
}
END_TEST

START_TEST( WriteVector_test )
{
  PetscErrorCode ierr;
  
  Vec v;
  VecCreate( PETSC_COMM_SELF, &v);
  VecSetSizes(v, 10, 10);
  VecSetType(v, VECSEQ );
  VecSet(v, one);
  
  ierr = WriteVector("testwrite.bin", v); CHKERRQ(ierr);
}
END_TEST

Suite* CreateSuite (void)
{
  Suite *s = suite_create ("Common Check");

  /* Core test case */
  TCase *tc_core = tcase_create("Core");
  tcase_add_checked_fixture( tc_core, setup, teardown);
//  tcase_add_test( tc_core, WriteVector_test );
//  tcase_add_test( tc_core, check_FloatAsIntArrayIndex );
//  tcase_add_test( tc_core, using_MKL_SVD );
//  tcase_add_test( tc_core, using_MKL_DOT );
//  tcase_add_test( tc_core, using_MKL_DGEMV );
//  tcase_add_test( tc_core, using_MKL_DROT  );
//  tcase_add_test( tc_core,  );  
//  tcase_add_test( tc_core, test_To1D );
  tcase_add_test( tc_core,  test_LeastSq );
//  tcase_add_test( tc_core,  test_PetscTemp );
//  tcase_add_test( tc_core,  test_WriteVec );
  suite_add_tcase( s, tc_core);

	TCase *tc_grid = tcase_create("Grid");
//	tcase_add_test( tc_grid, CreateGrid2D_test );
//	tcase_add_test( tc_grid, test_GridFunction );
  tcase_add_test( tc_grid, test_Grid3D );
  tcase_add_test( tc_grid, test_Bilinear3D );
	//  tcase_add_test( tc_grid,  );
	suite_add_tcase( s, tc_grid);
  return s;
}

int RunCheck()
{ 
  Suite *s = CreateSuite ();
  SRunner *sr = srunner_create (s);
  srunner_set_fork_status(sr, CK_NOFORK);
  srunner_run_all (sr, CK_NORMAL);
  srunner_free (sr);
  return 0;
}

void setup()
{
  PetscErrorCode ierr;
}

void teardown()
{
  PetscErrorCode ierr;
}

#undef __FUNCT__
#define __FUNCT__ "PetscMain"
PetscErrorCode PetscMain()
{
  PetscErrorCode ierr;
    
  PetscFunctionBegin;
  
  PetscFunctionReturn(0);
}
