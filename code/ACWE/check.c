#include "ActiveContoursWithoutEdges.h"
#include "MyCheck.h"

START_TEST( test_SingleStep )
{
  PetscErrorCode ierr;
  char *file_format = "/home/share/Images/Joanne/DIR.IGF1 Adaptation-10.spl/image.a.%d.gray";
  char *file_name;
  int file_name_len = 256;
  PetscMalloc( file_name_len * sizeof(char), &file_name);
  PetscSNPrintf(file_name,file_name_len,file_format, 10);
  
  
  Grid2D img;
  
  ierr = ReadCascade512File(file_name, &img); CHKERRQ(ierr);
  
  UserContext *uc;
  
  ierr = UserContextCreate( img, &uc); CHKERRQ(ierr);
  
  printf("%s\n", file_name );
  
  mark_point();
  ierr = SingleStep( uc, img); CHKERRQ(ierr);
  
}
END_TEST

Suite* CreateSuite (void)
{
  Suite *s = suite_create ("Common Check");

  /* Core test case */
  TCase *tc_core = tcase_create("Core");
  tcase_add_test( tc_core, test_SingleStep );
//  tcase_add_test( tc_core,  );  
  suite_add_tcase( s, tc_core);
  return s;
}

int RunCheck()
{ 
  Suite *s = CreateSuite ();
  SRunner *sr = srunner_create (s);
  srunner_set_fork_status(sr, CK_FORK);
  srunner_run_all (sr, CK_NORMAL);
  srunner_free (sr);
  return 0;
}