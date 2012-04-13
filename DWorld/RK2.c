#include "DWorld.h"

#undef __FUNCT__
#define __FUNCT__ "DWorldSimulate_RK2"
PetscErrorCode DWorldSimulate_RK2(DWorld w)
{
  FluidField fluid = w->fluid;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  for ( w->ti = 0; w->t < w->tend && w->ti < w->timax; ++w->ti) {
    // Reset rhs for next time step
    ierr = VecZeroEntries(fluid->rhs); CHKERRQ(ierr);
    ierr = DCellsArrayUpdateFluidFieldRHS(w->dcells, w->iim, fluid, w->t ); CHKERRQ(ierr);
    ierr = FluidFieldSolve( fluid ); CHKERRQ(ierr);
    ierr = FluidFieldMaxVelocityMag( fluid,&w->maxVel); CHKERRQ(ierr);

    // Calc dt
    w->dtcfl = w->CFL * fluid->dh.x / w->maxVel;
    w->dt = PetscMin(w->dtmax,w->dtcfl);
    if( w->dtframe > 0 && w->tiframe * w->dtframe <= w->t + w->dt ) {
      printf("\n%e <= %e\n", w->tiframe * w->dtframe, w->t + w->dt);
      w->dt = w->tiframe * w->dtframe - w->t;
      w->tiframe++;
    }

    if( w->dtframe < 0 ) {
      if( w->ti % w->writeInterval == 0 ) {
        ierr = DWorldWrite( w, w->ti); CHKERRQ(ierr);
      }
    } else {
      if( PetscAbs(w->t - (w->tiframe-1) * w->dtframe) < 1e-10 ) {
        ierr = DWorldWrite( w, w->tiframe); CHKERRQ(ierr);
      }
    }

    ierr = DCellsArrayAdvectRK2HalfStep( w->dcells, w->fluid, w->dt / 2. ); CHKERRQ(ierr);
    w->t = w->t + w->dt/2.;
    ierr = VecZeroEntries(fluid->rhs); CHKERRQ(ierr);
    ierr = DCellsArrayUpdateFluidFieldRHS(w->dcells, w->iim, fluid, w->t ); CHKERRQ(ierr);
    ierr = FluidFieldSolve( fluid ); CHKERRQ(ierr);
    ierr = DCellsArrayAdvectRK2FullStep( w->dcells, w->fluid, w->dt ); CHKERRQ(ierr);
    w->t = w->t + w->dt/2.;

    if( w->printStep ) {
      ierr = DWorldPrintStep(w); CHKERRQ(ierr);
    }

    if( w->ti == 0 ) {
      ierr = PetscLogStageRegister("Simulation loop", &w->stageSimLoop); CHKERRQ(ierr);
      ierr = PetscLogStagePush(w->stageSimLoop); CHKERRQ(ierr);
    }
  }
  ierr = PetscLogStagePop(); CHKERRQ(ierr);
  ierr = DWorldPrintStep(w); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
