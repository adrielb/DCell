#include "LevelSetMethod.h"

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc, char **args)
{
  PetscErrorCode  ierr;
  PetscFunctionBegin;
  ierr = DCellInit(); CHKERRQ(ierr);

  PetscReal radius = 1;
  PetscReal dx = 0.1;
  Coor dh = { dx, dx, dx};
  Coor center = { 0, 0, 0};
  
  LevelSet ls;
  ierr = LevelSetInitializeToSphere( dh, center, radius, &ls ); CHKERRQ(ierr);
  
  ierr = ArrayWrite(ls->band, "band", 0); CHKERRQ(ierr);
  ierr = LevelSetWriteIrregularNodeList( ls, 0); CHKERRQ(ierr);
  
  int len = ArrayLength(ls->irregularNodes);
  printf("len: %d\n", len);

  ierr = GridWrite( ls->phi, 0); CHKERRQ(ierr);
  
  ierr = DCellFinalize(); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
