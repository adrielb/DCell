#include "ImmersedInterfaceMethod.h"

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc, char **args)
{
  PetscErrorCode  ierr;
  ierr = DCellInit(); CHKERRQ(ierr);
  


  ierr = DCellFinalize(); CHKERRQ(ierr);
  return 0;
}
