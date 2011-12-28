#include "MicrofluidicSimulator.h"
#include "MyCheck.h"


START_TEST( viz_Matrices )
{
  PetscErrorCode ierr;
  
  
}
END_TEST

Suite* PressureIncrement_Suite(void)
{
  Suite *s = suite_create ("Pressure Increment suite");
  TCase *tc_core = tcase_create("Core");
  tcase_add_test( tc_core,  viz_Matrices );  
  //  tcase_add_test( tc_core,  );  
  tcase_set_timeout(tc_core, 0);
  suite_add_tcase( s, tc_core);

  return s;
}