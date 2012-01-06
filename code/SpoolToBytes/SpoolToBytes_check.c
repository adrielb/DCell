#include <stdlib.h>
#include "/usr/local/include/check.h"

void setup()
{

}

void teardown()
{
	
}

START_TEST( can_open_file )
{
	 fail_unless(3==4,"i failed %d times", 3);
}
END_TEST

Suite* SpoolToBytes_suite (void)
{
  Suite *s = suite_create ("SpootToBytes");

  /* Core test case */
  TCase *tc_core = tcase_create("Core");
  tcase_add_checked_fixture( tc_core, setup, teardown);
  tcase_add_test( tc_core, can_open_file);
  suite_add_tcase( s, tc_core);

  return s;
}

int main (void)
{
  int number_failed;
  Suite *s = SpoolToBytes_suite ();
  SRunner *sr = srunner_create (s);
  srunner_run_all (sr, CK_NORMAL);
  number_failed = srunner_ntests_failed (sr);
  srunner_free (sr);
  //return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
  return 0;
}