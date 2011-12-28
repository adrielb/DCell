#include "petsc.h"
#include "Main.h"
#include "MyCheck.h"


Suite* FluidField_Suite(void);

Suite* CreateSuite (void)
{
  Suite *s = suite_create ("Check");
  TCase *tc_core = tcase_create("Core");
  suite_add_tcase( s, tc_core);
  return s;
}

int RunCheck()
{ 
  Suite *s = CreateSuite ();
  SRunner *sr = srunner_create (s);
  srunner_add_suite(sr, FluidField_Suite() );
  srunner_set_fork_status(sr, CK_NOFORK);
  srunner_run_all (sr, CK_NORMAL);
  srunner_free (sr);
  return 0;
}

PetscErrorCode PetscMain()
{
  PetscFunctionReturn(0);
}