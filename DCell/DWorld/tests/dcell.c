#include "DWorld.h"

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc, char **args) {

  PetscErrorCode ierr;
  ierr = DCellInit(); CHKERRQ(ierr);

  DCell dcell;
  ierr = DCellCreate(&dcell); CHKERRQ(ierr);
  ierr = LevelSetCreate((Coor){1,1,1},(iCoor){0,0,0},(iCoor){10,10,0},&dcell->lsPlasmaMembrane); CHKERRQ(ierr);

  printf("dcell ID: %d\n", dcell->ID );

  ierr = dcell->Destroy(dcell); CHKERRQ(ierr);

	ierr = DCellFinalize(); CHKERRQ(ierr);
  return 0;
}
