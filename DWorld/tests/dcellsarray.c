#include "DWorld.h"

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc, char **args) {

  PetscErrorCode ierr;
  ierr = DCellInit(); CHKERRQ(ierr);

  DCellsArray dcells;
  ierr = DCellsArrayCreate(&dcells); CHKERRQ(ierr);

  DCell dcell;
  ierr = DCellCreate(&dcell); CHKERRQ(ierr);
  ierr = LevelSetInitializeToCircle( (Coor){1,1,1}, (Coor){0,0,0}, 10, &dcell->lsPlasmaMembrane); CHKERRQ(ierr);
  printf("len: %d\n", ArrayLength(dcells->dcells));
  ierr = DCellsArrayAdd(dcells, dcell); CHKERRQ(ierr);
  printf("len: %d\n", ArrayLength(dcells->dcells));
  ierr = DCellsArrayWrite(dcells,0); CHKERRQ(ierr);
  ierr = DCellsArrayDestroy(dcells); CHKERRQ(ierr);

	ierr = DCellFinalize(); CHKERRQ(ierr);
  return 0;
}
