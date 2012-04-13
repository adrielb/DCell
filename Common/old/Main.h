#ifndef MAIN_H
#define MAIN_H

#include "petsc.h"

int RunCheck();
PetscErrorCode PetscMain();

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc, char **args)
{
	PetscErrorCode  ierr;
	PetscTruth check;
	
	ierr = PetscInitialize(&argc, &args, (char *) 0, ""); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "Start\n"); CHKERRQ(ierr);
  
  PetscFunctionBegin;
  ierr = PetscOptionsHasName("", "-check", &check); CHKERRQ(ierr);
  
  if( check )
  	RunCheck();
  else
    PetscMain();
  
  ierr = PetscPrintf(PETSC_COMM_WORLD, "End\n"); CHKERRQ(ierr);
  ierr = PetscFinalize(); CHKERRQ(ierr);
 	PetscFunctionReturn(0);
}

#endif /* MAIN_H */