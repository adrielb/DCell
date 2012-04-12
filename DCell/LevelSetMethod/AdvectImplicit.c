#include "LevelSetMethod.h"
#include "LSM_private.h"

#undef __FUNCT__
#define __FUNCT__ "LevelSetAdvectImplicit"
PetscErrorCode LevelSetAdvectImplicit( LevelSet ls, Grid velgrid, PetscReal dt )
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = LevelSetAdvectSL(ls, velgrid, dt); CHKERRQ(ierr);
  ierr = LevelSet_MaxVelocity(ls, velgrid); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

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
  PetscReal *phi;
  PetscReal *psi;
  const PetscReal eps = 2.1;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = LevelSetAdvect(ls,ga,dt); CHKERRQ(ierr);
  ierr = GridGet(ls->phi,&phi); CHKERRQ(ierr);
//  ierr = GridGet(ls->tmp,&tmp); CHKERRQ(ierr);
//  ierr = VecZeroEntries(ls->tmp->v); CHKERRQ(ierr);
  ierr = GridGet(ls->psi->phi,&psi); CHKERRQ(ierr);
  if( ls->phi->is2D ) {
    PetscReal diracdelta;
    PetscReal **phi2D = (PetscReal**) phi;
    PetscReal **psi2D = (PetscReal**) psi;
    for (i = 0; i < len; ++i) {
      diracdelta = PetscAbs( phi2D[b[i].y][b[i].x] ) < eps ;
      g[i] = ( psi2D[b[i].y][b[i].x] - phi2D[b[i].y][b[i].x] ) * diracdelta;
//    tmp[b[i].y][b[i].x] = g[i];
    }
  } else {
    PetscReal diracdelta;
    PetscReal ***phi3D = (PetscReal***) phi;
    PetscReal ***psi3D = (PetscReal***) psi;
    for (i = 0; i < len; ++i) {
      diracdelta = PetscAbs( phi3D[b[i].z][b[i].y][b[i].x] ) < eps ;
      g[i] = ( psi3D[b[i].z][b[i].y][b[i].x] - phi3D[b[i].z][b[i].y][b[i].x] ) * diracdelta;
    }
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
  PetscReal *psi;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = GridGet(ls->psi->phi,&psi); CHKERRQ(ierr);
  if( ls->phi->is2D ) {
    PetscReal **psi2D = (PetscReal**)psi;
    for (i = 0; i < len; ++i) {
      psi2D[b[i].y][b[i].x] = psi2D[b[i].y][b[i].x] + lambda * dpsi[i];
    }
  } else {
    PetscReal ***psi3D = (PetscReal***)psi;
    for (i = 0; i < len; ++i) {
      psi3D[b[i].z][b[i].y][b[i].x] = psi3D[b[i].z][b[i].y][b[i].x] + lambda * dpsi[i];
    }
  }

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "LevelSetAdvectImplicitReinit"
PetscErrorCode LevelSetAdvectImplicitReinit( LevelSet ls, PetscReal dt )
{
  Coor dh = ls->phi->d;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = LevelSetUpdateIrregularNodeList( ls ); CHKERRQ(ierr);
  ls->CFLcount += dt * ls->maxVel / dh.x;
  ls->AdvectCount++;

  if( ls->AdvectCount > ls->AdvectThres || ls->CFLcount > ls->CFLthres ) {
    // TODO: specifically say which level set object is reinitializing (ls->ID)
    ierr = PetscInfo4(0,"CFLcount: %f:%f  AdvectCount: %d:%d\n", ls->CFLcount, ls->CFLthres, ls->AdvectCount,ls->AdvectThres); CHKERRQ(ierr);
    ls->CFLcount = 0;
    ls->AdvectCount = 0;
    ierr = LevelSetReinitialize(ls); CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}
