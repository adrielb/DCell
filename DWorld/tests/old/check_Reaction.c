#include "Reaction.h"
#include "MyCheck.h"

START_TEST( test_ComputeJacobian )
{
  PetscInt dof = 2;;
  PetscReal chem[2];
  Reaction rxn;
  PetscErrorCode ierr;
  ReactionCreate_Test(dof, &rxn);
  for (int i = 0; i < dof*dof; ++i) {
    fail_unless( rxn->jac[i] == i+1, "Jacobian evaluation error: %f", rxn->jac[i] );
  }
}
END_TEST

START_TEST( test_Reactions )
{
  PetscInt dof = 3, i;
  PetscReal chem[3]={1,2,3};
  Reaction rxn;
  PetscErrorCode ierr;
  
  ierr = ReactionCreate(dof, &rxn); CHKERRQ(ierr);
  mark_point();
  
  ReactionUpdateFunction(rxn, chem);
  for( i = 0; i < dof; ++i)
  {
    fail_unless(rxn->F[i] == 0, "Reaction not null: %f", rxn->F[i] );
    fail_unless(rxn->D[i] == 0, "Diffusion not zero");
  }

  ReactionUpdateJacobian(rxn, chem);
  for( i = 0; i < dof*dof; ++i)
  {
    fail_unless(rxn->jac[i] == 0, "Jacobian function not null");
  }
  
  ierr = ReactionDestroy(rxn); CHKERRQ(ierr);
}
END_TEST

Suite *Reaction_suite()
{
  Suite *s = suite_create ("Reaction Check");
  TCase *tc_core = tcase_create("Core");
  tcase_add_test( tc_core, test_Reactions );
  tcase_add_test( tc_core,  test_ComputeJacobian );
//  tcase_add_test( tc_core,  );
  tcase_set_timeout(tc_core, 0);
  suite_add_tcase( s, tc_core);
  return s;
}