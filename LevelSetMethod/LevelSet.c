#include "LevelSetMethod.h"
#include "LSM_private.h"

#undef __FUNCT__
#define __FUNCT__ "LevelSetCreate"
PetscErrorCode LevelSetCreate(Coor dh, iCoor pos, iCoor size, LevelSet *levelset)
{
  PetscReal bandWidth = 6;
  LevelSet ls;
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  ierr = PetscNew(struct _LevelSet, &ls); CHKERRQ(ierr);
  ierr = GridCreate(dh,pos,size,1,&ls->phi); CHKERRQ(ierr);
  ierr = GridSetName(ls->phi,"phi"); CHKERRQ(ierr);
  ierr = GridCreate(dh,pos,size,1,&ls->phi0); CHKERRQ(ierr);
  ierr = GridSetName(ls->phi0,"phi0"); CHKERRQ(ierr);
  ierr = GridCreate(dh,pos,size,1,&ls->tmp); CHKERRQ(ierr); // Use a global, shared buffer for all level sets
  ierr = GridSetName(ls->tmp,"tmp"); CHKERRQ(ierr);
  //TODO: set band width, init mem sizes for band and heap, counts and thresholds as petsc options
  ls->CFLthres = 3;
  ls->CFLcount = 0;
  ls->AdvectCount = 0;
  ls->AdvectThres = 1000;
  ierr = PetscOptionsGetReal(0,"-levelset_bandwidth",&bandWidth,0); CHKERRQ(ierr);
  ierr = LevelSetSetBandWidth(ls, bandWidth); CHKERRQ(ierr);
  // TODO: make a better estimate for initial sizes
  ierr = ArrayCreate( "irregularNodes", sizeof(IrregularNode), &ls->irregularNodes); CHKERRQ(ierr);
  ierr = ArrayCreate( "band", sizeof(iCoor), &ls->band); CHKERRQ(ierr);

  ls->Advect = LevelSetAdvectAndReinit;
  ierr = LevelSetRegisterEvents(); CHKERRQ(ierr);
  *levelset = ls;
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "LevelSetDestroy"
PetscErrorCode LevelSetDestroy(LevelSet ls)
{
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  if( ls->pls ) {
    ierr = ParticleLSDestroy(ls->pls); CHKERRQ(ierr);
  }
  ierr = GridDestroy(ls->phi); CHKERRQ(ierr);
  ierr = GridDestroy(ls->phi0); CHKERRQ(ierr);
  ierr = GridDestroy(ls->tmp); CHKERRQ(ierr);
  ierr = ArrayDestroy(ls->irregularNodes); CHKERRQ(ierr);
  ierr = ArrayDestroy(ls->band); CHKERRQ(ierr);
  ierr = PetscFree(ls); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "LevelSetSetBandWidth"
PetscErrorCode LevelSetSetBandWidth(LevelSet ls, PetscReal bandwidth)
{
  PetscFunctionBegin;
  ls->bandWidth = bandwidth;
  ls->PHI_INF = ls->bandWidth + 1.5; // bandWidth + sqrt(2) + eps
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "LevelSetDuplicate"
PetscErrorCode LevelSetDuplicate( LevelSet ls, LevelSet *copy)
{
  Grid phi = ls->phi;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = LevelSetCreate(phi->d,phi->p,phi->n,copy); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "LevelSetRegisterEvents"
static int EVENTS_registered = PETSC_FALSE;
PetscErrorCode LevelSetRegisterEvents(  )
{
  PetscErrorCode ierr;
  if( EVENTS_registered )
    PetscFunctionReturn(0);

  ierr = PetscLogEventRegister("LSWriteIrreg",0,&EVENT_LevelSetWriteIrregularNodeList); CHKERRQ(ierr);
  ierr = PetscLogEventRegister("LSUpdateIrreg",0,&EVENT_LevelSetUpdateIrregularNodeList); CHKERRQ(ierr);
  ierr = PetscLogEventRegister("LSReinit",0,&EVENT_LevelSetReinitialize); CHKERRQ(ierr);
  ierr = PetscLogEventRegister("LSGetVelGA",0,&EVENT_LevelSetGetVelocity); CHKERRQ(ierr);
  ierr = PetscLogEventRegister("LSAdvectSL",0,&EVENT_LevelSetAdvectSL); CHKERRQ(ierr);
  ierr = PetscLogEventRegister("LSAdvectParticle",0,&EVENT_ParticleLS_AdvectParticles); CHKERRQ(ierr);
  ierr = PetscLogEventRegister("LSAdjustParticle",0,&EVENT_ParticleLS_AdjustRadii); CHKERRQ(ierr);
  ierr = PetscLogEventRegister("LSErrorCorr",0,&EVENT_ParticleLS_ErrorCorrection); CHKERRQ(ierr);
  ierr = PetscLogEventRegister("LSReseedParticle",0,&EVENT_ParticleLS_ReseedParticles); CHKERRQ(ierr);

  EVENTS_registered = PETSC_TRUE;
  PetscFunctionReturn(0);
}
