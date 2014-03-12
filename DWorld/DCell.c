#include "DWorld.h"

static int DCELL_LOCAL_COUNT = 0;

PetscErrorCode DCellSetLocalCount( int start )
{
  DCELL_LOCAL_COUNT = start;
  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "DCellCreate"
PetscErrorCode DCellCreate( LevelSet lsPlasmaMembrane, DCell *cell )
{
  DCell c;
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  ierr = PetscNew( struct _DCell, &c); CHKERRQ(ierr);
  ierr = DCellSetup( lsPlasmaMembrane, c ); CHKERRQ(ierr);

  *cell = c;
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DCellSetup"
PetscErrorCode DCellSetup( LevelSet lsPlasmaMembrane, DCell cell )
{
  int size;
  int rank;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  cell->lsPlasmaMembrane = lsPlasmaMembrane;

  // Calculate unique ID for this cell
  ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size); CHKERRQ(ierr);
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank); CHKERRQ(ierr);
  cell->id = rank + size * DCELL_LOCAL_COUNT; // global ID of dcell
  DCELL_LOCAL_COUNT++;

  cell->type = 0;
  cell->Write = DCellWrite;
  cell->Destroy = DCellDestroy;
  cell->Advect = DCellAdvect;
  cell->AdvectRK2HalfStep = DCellAdvectRK2HalfStep;
  cell->AdvectRK2FullStep = DCellAdvectRK2FullStep;
  cell->UpdateFluidFieldRHS = DCellUpdateFluidFieldRHS;
  cell->AdvectImplicitInit = DCellAdvectImplicitInit;
  cell->AdvectImplicitRHS = DCellAdvectImplicitRHS;
  cell->AdvectImplicitUpdate = DCellAdvectImplicitUpdate;
  cell->AdvectImplicitReinit = DCellAdvectImplicitReinit;

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DCellDestroy"
PetscErrorCode DCellDestroy( DCell c )
{
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  ierr = LevelSetDestroy(c->lsPlasmaMembrane); CHKERRQ(ierr);
  ierr = PetscFree(c); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DCellWrite"
PetscErrorCode DCellWrite( DCell dcell, int ti )
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = GridWrite(dcell->lsPlasmaMembrane->phi, ti); CHKERRQ(ierr);
  ierr = LevelSetWriteIrregularNodeList(dcell->lsPlasmaMembrane, ti); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DCellAdvect"
PetscErrorCode DCellAdvect( DCell dcell, int ga, PetscReal dt )
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = LevelSetAdvect(dcell->lsPlasmaMembrane, ga, dt); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DCellAdvectRK2HalfStep"
PetscErrorCode DCellAdvectRK2HalfStep( DCell dcell, int ga, PetscReal dt )
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  dcell->lsPlasmaMembrane->Advect = LevelSetAdvectPLSRK2HalfStep;
  ierr = LevelSetAdvect(dcell->lsPlasmaMembrane, ga, dt); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DCellAdvectRK2FullStep"
PetscErrorCode DCellAdvectRK2FullStep( DCell dcell, int ga, PetscReal dt )
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  dcell->lsPlasmaMembrane->Advect = LevelSetAdvectPLSRK2FullStep;
  ierr = LevelSetAdvect(dcell->lsPlasmaMembrane, ga, dt); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DCellUpdateFluidFieldRHS"
PetscErrorCode DCellUpdateFluidFieldRHS( DCell this, IIM iim, int ga, PetscReal t )
{
  PetscFunctionBegin;
  
  //TODO: implement DCellUpdateFluidFieldRHS
  if (t && this && iim && ga) {
  }
  SETERRQ(PETSC_COMM_SELF,0,"Not Implemented");

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DCellAdvectImplicitInit"
PetscErrorCode DCellAdvectImplicitInit( DCell this, PetscInt *n )
{
  PetscErrorCode ierr;
  PetscFunctionBegin;
  ierr = LevelSetAdvectImplicitInit( this->lsPlasmaMembrane, n); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DCellAdvectImplicitRHS"
PetscErrorCode DCellAdvectImplicitRHS( DCell this, int ga, PetscReal dt, PetscReal *g )
{
  PetscErrorCode ierr;
  PetscFunctionBegin;
  ierr = LevelSetAdvectImplicitRHS( this->lsPlasmaMembrane, ga, dt, g ); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DCellAdvectImplicitUpdate"
PetscErrorCode DCellAdvectImplicitUpdate( DCell this, PetscReal lambda, PetscReal *dpsi )
{
  PetscErrorCode ierr;
  PetscFunctionBegin;
  ierr = LevelSetAdvectImplicitUpdate( this->lsPlasmaMembrane, lambda, dpsi ); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DCellAdvectImplicitReinit"
PetscErrorCode DCellAdvectImplicitReinit( DCell this, PetscReal dt )
{
  PetscErrorCode ierr;
  PetscFunctionBegin;
  ierr = LevelSetAdvectImplicitReinit( this->lsPlasmaMembrane, dt ); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
