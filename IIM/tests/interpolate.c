#include "ImmersedInterfaceMethod.h"

inline void IIMVelocityCorrection( Coor X, IrregularNode *n, Coor *vel );

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc, char **args)
{
  PetscErrorCode  ierr;
  ierr = DCellInit(); CHKERRQ(ierr);

  IrregularNode n = {
      .shift = CELL_CENTER,
//      .shift = U_FACE,
      .axis  = 0,
      .X     = {0.8,0,0},
      .pos   = {1,0,0},
      .uj    = 10
  };

  PetscViewer velfield;
  char filename[] = "velfield.dat";
  ierr = PetscViewerASCIIOpen(PETSC_COMM_SELF,filename,&velfield); CHKERRQ(ierr);

  Coor X;
  for (X.y = -1.5; X.y < 1.5; X.y += 0.02 ) {
    for (X.x = -1.5; X.x < 1.5; X.x += 0.02 ) {
      Coor vel = {0,0,0};
      PetscReal *v = &vel.x;
      IIMVelocityCorrection( X, &n, &vel);
      ierr = PetscViewerASCIIPrintf(velfield,"%f ", v[0] ); CHKERRQ(ierr);
    }
  }

  ierr = PetscViewerDestroy(&velfield); CHKERRQ(ierr);
  ierr = DCellFinalize(); CHKERRQ(ierr);
  return 0;
}
