#include "LevelSetMethod.h"
#include "MyCheck.h"
#include "Main.h"

START_TEST( viz_SphereInit )
{
  PetscErrorCode ierr;
  int d1=128, d2=128, d3=16;
  LevelSet3D ls;
  LevelSet3DCreate("SphereInit", d1, d2, d3, &ls);
  LevelSetInitializeToSphere( ls, PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE );
  WriteVecToDisk(ls->g3d->wv);
}
END_TEST

START_TEST( viz_StarInit )
{
  PetscErrorCode ierr;
  int d1=300, d2=300;
  LevelSet2D ls;
  CreateLevelSet2D(d1,d2,&ls);
  LevelSetInitializeToStar(ls,PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE);
  WriteVector("StarInit",ls->g2d->v);
}
END_TEST

START_TEST( test_OrthoProjFirstOrder2D )
{
  double p[9] = {0, 1, 0,
                 1,-1, 1,
                 0, 1, 0};
  double u, v, mindist;
  int i=1, j=1;
  PetscErrorCode ierr;
  
  Grid2D g;
  CreateGrid2D(3,3,&g);
  for( int i = 0; i<9; i++)
    g->v1[i] = p[i];
    
  fail("problem with args:");
//  ierr = OrthoProjFirstOrder2D(g->v2, i, j, &u, &v, &mindist); CHKERRQ(ierr);
  
//  printf("%f, %f\t%f\n", u, v, mindist);
}
END_TEST

//TODO reinit needs to set "regular" grid points to INF

START_TEST( print_ReinitializeLevelSet )
{
  LevelSet2D ls;
  
  int d1 = 8, d2 = 8, dd=64;
  double d[64] = {
    1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1,-1,-1, 1, 1, 1,
    1, 1, 1,-1,-1, 1, 1, 1,
    1, 1, 1,-1,-1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1
  };
  PetscErrorCode ierr;
  
  CreateLevelSet2D(d1,d2,&ls);
  double *p = ls->g2d->v1, **phi = ls->g2d->v2;
  
  for( int i = 0; i < dd; i++ )
    p[i] = d[i];
  
  UpdateIrregularNodeList(ls);
  
  IrregularNode n;
  for( int i = 0; i < ls->irregularNodes->len; i++ )
  {
    n = g_array_index(ls->irregularNodes, IrregularNode, i);
    printf("%f\t%d,%d,%d\n\t%f,%f,%f\n",n.sign,n.x, n.y, n.z,n.ox,n.oy,n.oz);
  }  

  ReinitializeLevelSet( ls );
  
  PrintPhi(d1,d2,phi);
  
  DestroyLevelSet2D(ls);
}
END_TEST

START_TEST( print_IndexBCNodes )
{
  LevelSet2D ls;
  
  int d1 = 8, d2 = 8, dd=64;
  double d[64] = {
    1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1,-1, 1, 1, 1, 1,
    1, 1, 1,-1,-1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1
  };
  PetscErrorCode ierr;
  
  CreateLevelSet2D(d1,d2,&ls);
  double *p = ls->g2d->v1, **phi = ls->g2d->v2;
  
  for( int i = 0; i < dd; i++ )
    p[i] = d[i];
  
  UpdateIrregularNodeList( ls );
  
  IrregularNode n;
  for( int i = 0; i< ls->irregularNodes->len; i++)
  {
    n = g_array_index(ls->irregularNodes, IrregularNode, i);
    printf( "P: %d, %d\n", n.x, n.y);
  }
  
  BCFromOrthoProj( ls );
  
  PrintPhi(d1, d2, phi);
}
END_TEST

START_TEST( one_dim_to_two_dim )
{
  PetscErrorCode ierr;
  
  int d1 = 4, d2 = 6, d, c1, c2;
  
  for( int i = 0; i < d1; i++ )
  {
    for( int j = 0; j < d2; j++ )
    {
      d = i + j * d1;
      c2 = d / d1;
      c1 = d - d1 * c2;
      fail_unless(c1 == i && c2 == j, "c1 = %d\ti=%d\tc2 = %d\tj=%d", c1, i, c2, j);
    } 
  }
}
END_TEST

START_TEST( test_GArray_intArray )
{
  PetscErrorCode ierr;
  
  GArray *garray = g_array_new( FALSE, FALSE, 2 * sizeof(int) );
  {
    Coor2D c;
    for( int i = 0; i < 10; i++ )
    {
      c.x = i;
      c.y = -1 * i;
      g_array_append_val(garray, c);
    }
  }
  
  Coor2D c;
  
  for( int i = 0; i < 10; i++)
  {
    c = g_array_index(garray, Coor2D, i);
    fail_unless( c.x == i, "c[0] == %d, i: %d", c.x, i);
    fail_unless( c.y == -1 * i, "c[0] == %d, i: -%d", c.y, i);
  }
  
}
END_TEST

START_TEST( test_GArray )
{
  int i;
  PetscErrorCode ierr;
  
  GArray *garray = g_array_new(FALSE, FALSE, sizeof(int) );
  
  // Append ints 1-10 into garray
  for( i = 0; i < 10; i++)
    g_array_append_val(garray, i);

  // Check g_array_index gets ints 1-10
  for (i = 0; i < 10; i++)
    fail_unless( g_array_index (garray, int, i) == i,
              "ERROR: got %d instead of %d",
               g_array_index (garray, int, i), i);
  
  // Check using private data member 
  int *data = (int*)garray->data;
  for (i = 0; i < 10; i++)
    fail_unless( data[i] == i,
                 "ERROR: got %d instead of %d",
                 data[i], i);
  g_array_free( garray, TRUE);
}
END_TEST

Suite* ZD_suite (void);
Suite* FMM_suite (void)
{
  Suite *s = suite_create ("Fast Marching Method Check");
  TCase *tc_core = tcase_create("Core");
  tcase_add_checked_fixture( tc_core, setup, teardown);
  tcase_add_test( tc_core, test_GArray);
  tcase_add_test( tc_core, test_GArray_intArray );
  tcase_add_test( tc_core, one_dim_to_two_dim );
//  tcase_add_test( tc_core, print_IndexBCNodes );
//  tcase_add_test( tc_core, print_ReinitializeLevelSet );
  tcase_add_test( tc_core,  test_OrthoProjFirstOrder2D );
  tcase_add_test( tc_core,  viz_StarInit );
  tcase_add_test( tc_core,  viz_SphereInit );
//  tcase_add_test( tc_core,  );
  suite_add_tcase( s, tc_core);
  return s;
}

int RunCheck()
{ 
  Suite *s = FMM_suite ();
//  SRunner *sr = srunner_create (ZD_suite());
  SRunner *sr = srunner_create (s);
//  srunner_add_suite(sr, s );
  srunner_set_fork_status(sr, CK_FORK);
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

PetscErrorCode PetscMain() 
{
  PetscFunctionReturn(0);
}