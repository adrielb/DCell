#include "LevelSetMethod.h"

#undef __FUNCT__
#define __FUNCT__ "LevelSetAdvectSLRK2HalfStep"
PetscErrorCode LevelSetAdvectSLRK2HalfStep( LevelSet ls, Grid velgrid, PetscReal dt )
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = GridCopy(ls->phi,ls->phi0); CHKERRQ(ierr);
  ierr = LevelSetAdvectSL(ls, velgrid, dt ); CHKERRQ(ierr);
  ierr = LevelSetUpdateIrregularNodeList( ls ); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "LevelSetAdvectSLRK2FullStep"
PetscErrorCode LevelSetAdvectSLRK2FullStep( LevelSet ls, Grid velgrid, PetscReal dt )
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = LevelSetAdvectSL(ls, velgrid, dt ); CHKERRQ(ierr);
  ierr = LevelSetUpdateIrregularNodeList( ls ); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

// Particle Level Set Advection

