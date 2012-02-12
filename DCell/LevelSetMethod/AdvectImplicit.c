#include "LevelSetMethod.h"
#include "LSM_private.h"

#undef __FUNCT__
#define __FUNCT__ "LevelSetAdvectImplicitInit"
PetscErrorCode LevelSetAdvectImplicitInit( LevelSet ls, PetscInt *n )
{
  int len = ArrayLength(ls->band);
  PetscErrorCode ierr;

  PetscFunctionBegin;

  *n = *n + len;

  if( ls->psi == NULL ) {
    ierr = LevelSetDuplicate(ls, &ls->psi); CHKERRQ(ierr);
    ierr = GridSetName(ls->psi->phi,"psi"); CHKERRQ(ierr);
  }
  ierr = GridCopy(ls->phi,ls->psi->phi); CHKERRQ(ierr);
  ierr = ArrayCopy(ls->band,ls->psi->band); CHKERRQ(ierr);
  ierr = GridCopy(ls->phi,ls->phi0); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "LevelSetAdvectImplicitRHS"
PetscErrorCode LevelSetAdvectImplicitRHS( LevelSet ls, int ga, PetscReal dt, PetscReal *g )
{
  int i;
  iCoor *b = ArrayGetData(ls->band);
  int len = ArrayLength(ls->band);
  PetscReal **phi;
  PetscReal **psi;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = LevelSetAdvect(ls,ga,dt); CHKERRQ(ierr);
  ierr = GridGet(ls->phi,&phi); CHKERRQ(ierr);
  ierr = GridGet(ls->psi->phi,&psi); CHKERRQ(ierr);
  for (i = 0; i < len; ++i) {
    g[i] = psi[b[i].y][b[i].x] - phi[b[i].y][b[i].x];
  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "LevelSetAdvectImplicitUpdate"
PetscErrorCode LevelSetAdvectImplicitUpdate( LevelSet ls, PetscReal lambda, PetscReal *dpsi )
{
  int i;
  iCoor *b = ArrayGetData(ls->band);
  int len = ArrayLength(ls->band);
  PetscReal **psi;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = GridGet(ls->psi->phi,&psi); CHKERRQ(ierr);
  for (i = 0; i < len; ++i) {
    psi[b[i].y][b[i].x] = psi[b[i].y][b[i].x] + lambda * dpsi[i];
  }
  PetscFunctionReturn(0);
}
