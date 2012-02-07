#include "LevelSetMethod.h"
#include "LSM_private.h"

#undef __FUNCT__
#define __FUNCT__ "LevelSetAdvectImplicitRHS"
PetscErrorCode LevelSetAdvectImplicitRHS( LevelSet ls, Grid velgrid, PetscReal dt, PetscReal *p, PetscReal *g )
{
  int i;
  iCoor *b = ArrayGetData(ls->band);
  int len = ArrayLength(ls->band);
  PetscReal **phi;
  PetscReal **psi;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = LevelSetAdvectSL(ls,velgrid,dt); CHKERRQ(ierr);
  ierr = GridGet(ls->phi,&phi); CHKERRQ(ierr);
  ierr = GridGet(ls->psi,&psi); CHKERRQ(ierr);
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
  ierr = GridGet(ls->phi,&psi); CHKERRQ(ierr);
  for (i = 0; i < len; ++i) {
    psi[b[i].y][b[i].x] = psi[b[i].y][b[i].x] + lambda * dpsi[i];
  }
  PetscFunctionReturn(0);
}
