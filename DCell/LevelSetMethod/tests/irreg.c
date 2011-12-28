#include "LevelSetMethod.h"

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc, char **args)
{
  PetscErrorCode  ierr;
  PetscFunctionBegin;
  ierr = DCellInit(); CHKERRQ(ierr);

  PetscReal dx = 1;
  PetscReal radius = 10;
  Coor center = {12,12,0};
  Coor dh = {dx,dx,0};

  LevelSet ls;
//  ierr = LevelSetInitializeToCircle( dh, (Coor){0,0,0}, radius, &ls); CHKERRQ(ierr);
  ierr = LevelSetInitializeToStar2D(dh,center,radius,0.5*radius, 5,&ls); CHKERRQ(ierr);
  ierr = GridWrite(ls->phi,0); CHKERRQ(ierr);
  ierr = LevelSetWriteIrregularNodeList(ls->irregularNodes,0); CHKERRQ(ierr);
  ierr = LevelSetDestroy(ls); CHKERRQ(ierr);
  ierr = DCellFinalize(); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
