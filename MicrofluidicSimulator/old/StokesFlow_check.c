#include "MicrofluidicSimulator.h"
#include "MyCheck.h"

char *checkParams = "check.params";

START_TEST( DiffusiveAdvective_test )
{
  int i, j;
  PetscReal val[5];
  PetscReal DIFFUSION=1, DELTA_T=1, DELTA_X=1, ut, vt, u=-1, v=1;
  PetscInt numNei = 4;
  
  // Diffusive part
  for( j = 0; j < 4; j++)
    val[j] = (-1 * DIFFUSION * DELTA_T) / DELTA_X;
  val[4] = 1. + ( numNei * DIFFUSION * DELTA_T ) / ( DELTA_X * DELTA_X );
  
  // Advective part (upwinding) 
  ut = u * DELTA_T;
  vt = v * DELTA_T;
  if( u < 0. ) { val[1] += ut; val[4] -= ut * DELTA_X; }
  if( u > 0. ) { val[2] -= ut; val[4] += ut * DELTA_X; }
  if( v < 0. ) { val[0] += vt; val[4] -= vt * DELTA_X; }
  if( v > 0. ) { val[3] -= vt; val[4] += vt * DELTA_X; }
  
  printf("\nDiffusiveAdvective_test:\t");
  for( i = 0; i < 5; i++ )
    printf("%f\t", val[i] );
  printf("\n"); 
}
END_TEST

START_TEST( IndexFreeNodes_test )
{
  PetscErrorCode ierr;
  int i;
  UserContext uc_stack;
  UserContext *uc = &uc_stack;
  
  ierr = InterpretOptions(uc); CHKERRQ(ierr);
  ierr = ReadFile(uc->imageFileName, uc->n, &uc->filedata); CHKERRQ(ierr);
  printf("\nn:%d\n",uc->n);
  ierr = IndexFreeNodes(uc); CHKERRQ(ierr);
  printf("\n\n");
  for( i = 0; i < uc->n; i++)
    if( i%uc->numcols == 0 )
      printf("\n%d", (uc->filedata[i]) );
    else
      printf("\t%d", (uc->filedata[i]) );   
  
  printf("\n\n");   
  
  for( i = 0; i < uc->n; i++)
    if( i%uc->numcols == 0 )
      printf("\n%d",uc->imageToNode[i]);
    else
      printf("\t%d",uc->imageToNode[i]);
  
  printf("\n\n");
  
  for( i = 0; i < uc->numBC; i++)
    printf("BC %d: %d\n", i, uc->bcNodes[i].nodeIndex );
  
  printf("\n\n");
}
END_TEST

Suite* StokesFlow_suite (void)
{
  Suite *s = suite_create ("Stokes Flow Check");

  /* Core test case */
  TCase *tc_core = tcase_create("Core");
  tcase_add_checked_fixture( tc_core, setup, teardown);
  
  
//  tcase_add_test( tc_core, IndexFreeNodes_test );
  
//  tcase_add_test( tc_core, DiffusiveAdvective_test );
  
//  tcase_add_test( tc_core,  );  
  tcase_set_timeout(tc_core, 0);
  suite_add_tcase( s, tc_core);

  return s;
}

Suite* PressureIncrement_Suite(void);

int RunCheck()
{ 
  Suite *s = StokesFlow_suite ();
  SRunner *sr = srunner_create (s);
  srunner_add_suite(sr, PressureIncrement_Suite());
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