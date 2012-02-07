#include "DWorld.h"

PetscLogEvent EVENT_DWorldWrite;

#undef __FUNCT__
#define __FUNCT__ "DWorldCreate"
PetscErrorCode DWorldCreate( FluidField fluid, DWorld *world )
{
  struct _DWorld *w;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscNew(struct _DWorld, &w); CHKERRQ(ierr);
  
  char filename[PETSC_MAX_PATH_LEN];
  char wd[PETSC_MAX_PATH_LEN];
  ierr = DCellGetWorkingDirectory(wd); CHKERRQ(ierr);
  ierr = PetscSNPrintf(filename,PETSC_MAX_PATH_LEN,"%s/temporal.dat",wd); CHKERRQ(ierr);
  ierr = PetscViewerASCIIOpen(PETSC_COMM_SELF,filename,&w->temporalfile); CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(w->temporalfile,"iter t dt dtcfl tiframe CFL\n"); CHKERRQ(ierr);
  
  ierr = DCellsArrayCreate(&w->dcells); CHKERRQ(ierr);
  ierr = IIMCreate(!fluid->is3D,&fluid->mu,1.0,32,fluid->dh,&w->iim); CHKERRQ(ierr);

  w->fluid = fluid;
  w->timax = 1;
  w->tend = 1;
  w->CFL = 0.9;
  w->dt = 1;
  w->writeInterval = 1;
  w->dtframe = -1;
  w->tiframe = 0;
  w->printStep = PETSC_TRUE;
  w->Simulate = DWorldSimulate_BFGS;

  //BFGS
  int N = 100;
  MPI_Comm comm = PETSC_COMM_SELF;
  ierr = MatCreateSeqDense(comm,N,N,PETSC_NULL,&w->jac); CHKERRQ(ierr);
  ierr = VecCreateSeq(comm,N,&w->s); CHKERRQ(ierr);
  ierr = VecDuplicate(w->s,&w->y); CHKERRQ(ierr);
  ierr = VecDuplicate(w->s,&w->x0); CHKERRQ(ierr);
  ierr = VecDuplicate(w->s,&w->x1); CHKERRQ(ierr);

  ierr = PetscLogEventRegister("DWorldWrite", 0, &EVENT_DWorldWrite); CHKERRQ(ierr);
  ierr = PetscInfo(0, "Created DWorld\n"); CHKERRQ(ierr);

  *world = w;
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DWorldDestroy"
PetscErrorCode DWorldDestroy(DWorld world)
{
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  ierr = PetscViewerDestroy(&world->temporalfile); CHKERRQ(ierr);
  ierr = DCellsArrayDestroy(world->dcells); CHKERRQ(ierr);
  ierr = IIMDestroy(world->iim); CHKERRQ(ierr);
  ierr = FluidFieldDestroy(world->fluid); CHKERRQ(ierr);
  ierr = PetscFree(world); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DWorldAddDCell"
PetscErrorCode DWorldAddDCell( DWorld world, void *dcell )
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = DCellsArrayAdd( world->dcells, dcell ); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DWorldWrite"
PetscErrorCode DWorldWrite( DWorld world, int ti )
{
  char line[128];
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscLogEventBegin(EVENT_DWorldWrite,0,0,0,0); CHKERRQ(ierr);
  ierr = FluidFieldWrite( world->fluid,  ti); CHKERRQ(ierr);
  ierr = DCellsArrayWrite(world->dcells, ti); CHKERRQ(ierr);
  sprintf(line,"%d %e %e %e %d %e\n",
      world->ti,world->t,world->dt,world->dtcfl,world->tiframe, world->dt * world->maxVel / world->fluid->dh.x);
  ierr = PetscInfo1(0,"%s",line); CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(world->temporalfile,"%s",line); CHKERRQ(ierr);
  ierr = PetscLogEventEnd(EVENT_DWorldWrite,0,0,0,0); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "DWorldSimulate"
PetscErrorCode DWorldSimulate(DWorld w)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = w->Simulate( w ); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#include "petscblaslapack.h"
#undef __FUNCT__
#define __FUNCT__ "DWorldSimulate_BFGS"
PetscErrorCode DWorldSimulate_BFGS(DWorld w)
{
  int i, j, k;
  PetscReal zero = 0, one = 1, negone = -1;
  PetscReal sTy, yTBy, a1;
  PetscBLASInt n = 1;
  PetscReal *m1, *m2, *ssT, *BysT, *syTB, *B; // matrices
  PetscReal *s, *y, *d, *x0, *x1, *g0, *g1, *By; // vectors
  int inc = 1;
  const char nop = 'N';
  int lda;
  int maxIter = 10;
  PetscReal tol = 10;
  PetscReal lambda = 1;
  FluidField fluid = w->fluid;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  // set [B] = [I]
  // set phi0 = phi
  size_t vecsize = sizeof(double)*n;
  ierr = PetscMalloc(vecsize,&s); CHKERRQ(ierr);
  ierr = PetscMalloc(vecsize,&y); CHKERRQ(ierr);
  ierr = PetscMalloc(vecsize,&d); CHKERRQ(ierr);
  ierr = PetscMalloc(vecsize,&x0); CHKERRQ(ierr);
  ierr = PetscMalloc(vecsize,&x1); CHKERRQ(ierr);
  ierr = PetscMalloc(vecsize,&g0); CHKERRQ(ierr);
  ierr = PetscMalloc(vecsize,&g1); CHKERRQ(ierr);
  ierr = PetscMalloc(vecsize,&By); CHKERRQ(ierr);
  size_t matsize = sizeof(double)*n*n;
  ierr = PetscMalloc(matsize,&m1); CHKERRQ(ierr);
  ierr = PetscMalloc(matsize,&m2); CHKERRQ(ierr);
  ierr = PetscMalloc(matsize,&ssT); CHKERRQ(ierr);
  ierr = PetscMalloc(matsize,&BysT); CHKERRQ(ierr);
  ierr = PetscMalloc(matsize,&syTB); CHKERRQ(ierr);
  ierr = PetscMalloc(matsize,&B); CHKERRQ(ierr);

  for (i = 0; i < maxIter; ++i) {
    // Update RHS
    // Evaluate g0 = g( X = x0 )
    // g(X) = X - Xn - dt*U(X)
    ierr = VecZeroEntries(fluid->rhs); CHKERRQ(ierr);
    ierr = DCellsArrayUpdateFluidFieldRHS(w->dcells, w->iim, fluid, w->t ); CHKERRQ(ierr);
    ierr = FluidFieldSolve( fluid ); CHKERRQ(ierr);
    ierr = DCellsArrayAdvectImplicit( w->dcells, w->fluid, w->dt, x0, g0 ); CHKERRQ(ierr);

    // Solve B.d = -g0
    // Mult  d = -B.g0
    BLASgemv_(&nop,&n,&n,&negone,B,&lda,g0,&inc,&zero,d,&inc);

    // Update position
    // x1 = x0 + l d
    for (j = 0; j < n; ++j) {
      x1[j] = x0[j] + lambda*d[j];
    }

    // Check convergence
    // | x1 - x0 | < tol
    //  | g(x0) |  < tol
    PetscReal dx, norm;
    norm = 0;
    for (k = 0; k < n; ++k) {
      dx = x1[k] - x0[k];
      norm += dx*dx;
    }
    if( norm < tol )
      break;

    // Evaluate g1 = g( X = x1 )
    // g(X) = X - Xn - dt*U(X)
    ierr = VecZeroEntries(fluid->rhs); CHKERRQ(ierr);
    ierr = DCellsArrayUpdateFluidFieldRHS(w->dcells, w->iim, fluid, w->t ); CHKERRQ(ierr);
    ierr = FluidFieldSolve( fluid ); CHKERRQ(ierr);
    ierr = DCellsArrayAdvectImplicit( w->dcells, w->fluid, w->dt, x1, g1 ); CHKERRQ(ierr);

    // Update Hessian
    // s = x1 - x0
    // y = g(x1) - g(x0)
    // B = B + (sT.y + yT.B.y)(s.sT) / (sT.y)^2 - ( B.y.sT + s.yT.B ) / (sT.y)
    for (j = 0; j < n; ++j) {
      s[j] = lambda*d[j];
      y[j] = g1[j] - g0[j];
    }

    // [ (sT.y + yT.B.y) / (sT.y)^2 ] (s.sT)
    BLASgemv_(&nop,&n,&n,&one,B,&lda,y,&inc,&zero,By,&inc); // By = B.y
    yTBy = BLASdot_(&n,y,&inc,By,&inc);  // yT.B.y = y.By
    sTy  = BLASdot_(&n,s,&inc,y,&inc); // sTy = sT.y
    a1 = ( sTy + yTBy ) / (sTy*sTy);
    dger_(&n,&n,&a1,s,&inc,s,&inc,ssT,&lda); // ssT = a1 * s.sT
    // ( B.y.sT + s.yT.B ) / (sT.y)
    dger_(&n,&n,&one,y,&inc,s,&inc,m1,&lda); // m1 = y.sT
    dger_(&n,&n,&one,s,&inc,y,&inc,m2,&lda); // m2 = s.yT
    dgemm_(&nop,&nop,&n,&n,&n,&one,B,&lda,m1,&lda,&zero,BysT,&lda); // BysT = B.m1
    dgemm_(&nop,&nop,&n,&n,&n,&one,m2,&lda,B,&lda,&zero,syTB,&lda); // syTB = m2.B

    for (j = 0; j < n*n; ++j) {
      B[i] = B[i] + ssT[i] - (BysT[i] + syTB[i])/sTy;
    }
  }
  ierr = PetscFree(s); CHKERRQ(ierr);
  ierr = PetscFree(y); CHKERRQ(ierr);
  ierr = PetscFree(d); CHKERRQ(ierr);
  ierr = PetscFree(x0); CHKERRQ(ierr);
  ierr = PetscFree(x1); CHKERRQ(ierr);
  ierr = PetscFree(g0); CHKERRQ(ierr);
  ierr = PetscFree(g1); CHKERRQ(ierr);
  ierr = PetscFree(By); CHKERRQ(ierr);
  ierr = PetscFree(m1); CHKERRQ(ierr);
  ierr = PetscFree(m2); CHKERRQ(ierr);
  ierr = PetscFree(ssT); CHKERRQ(ierr);
  ierr = PetscFree(BysT); CHKERRQ(ierr);
  ierr = PetscFree(syTB); CHKERRQ(ierr);
  ierr = PetscFree(B); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

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

#undef __FUNCT__
#define __FUNCT__ "DWorldSimulate_Euler"
PetscErrorCode DWorldSimulate_Euler(DWorld w)
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

    w->dtcfl = w->CFL * fluid->dh.x / w->maxVel;
    w->dt = PetscMin(w->dtmax,w->dtcfl);
    if( w->dtframe < 0 ) {
      if( w->ti % w->writeInterval == 0 ) {
        ierr = DWorldWrite( w, w->ti); CHKERRQ(ierr);
      }
    } else {
      if( PetscAbs(w->t - (w->tiframe-1) * w->dtframe) < 1e-10 ) {
        ierr = DWorldWrite( w, w->tiframe); CHKERRQ(ierr);
      }
      if( w->tiframe * w->dtframe <= w->t + w->dt ) {
        printf("\n%e <= %e\n", w->tiframe * w->dtframe, w->t + w->dt);
        w->dt = w->tiframe * w->dtframe - w->t;
        w->tiframe++;
      }
    }

    if( w->printStep ) {
      ierr = DWorldPrintStep(w); CHKERRQ(ierr);
    }

    w->t = w->t + w->dt;

    ierr = DCellsArrayAdvect( w->dcells, w->fluid, w->dt ); CHKERRQ(ierr);

    if( w->ti == 0 ) {
      ierr = PetscLogStageRegister("Simulation loop", &w->stageSimLoop); CHKERRQ(ierr);
      ierr = PetscLogStagePush(w->stageSimLoop); CHKERRQ(ierr);
    }
  }
  ierr = PetscLogStagePop(); CHKERRQ(ierr);
  ierr = DWorldPrintStep(w); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode DWorldPrintStep( DWorld w )
{
  PetscErrorCode ierr;
  if( w->dtframe > 0 ) {
    ierr = PetscPrintf(PETSC_COMM_WORLD, "frame:\t %8d \n", w->tiframe-1);
  }
  ierr = PetscPrintf(PETSC_COMM_WORLD, "iter: \t %8d \t %8d \t %3.2f%%\n", w->ti, w->timax, (100.*w->ti)/w->timax);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "time: \t %.3e \t %.3e \t %3.2f%%\n", w->t, w->tend, (100.*w->t)/w->tend);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "dt:   \t %.3e \t %.3e \n", w->dt, w->dtcfl);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "CFL:  \t %f \t %f \n", w->dt*w->maxVel/w->fluid->dh.x, w->CFL);
  return 0;
}

PetscErrorCode DWorldSetPrintStep( DWorld w, PetscBool printStep )
{
  w->printStep = printStep;
  return 0;
}

PetscErrorCode DWorldSetWriteInterval( DWorld w, PetscInt interval )
{
  w->writeInterval = interval;
  return 0;
}

PetscErrorCode DWorldSetFrameInterval( DWorld w, PetscReal dtframe )
{
  w->dtframe = dtframe;
  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "DWorldSetFromOptions"
PetscErrorCode DWorldSetFromOptions( DWorld w )
{
  MPI_Comm comm = w->fluid->comm;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscOptionsGetReal(0,"-CFL",&w->CFL,0); CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(0,"-dtmax",&w->dtmax,0); CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(0,"-dtframe",&w->dtframe,0); CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(0,"-tend",&w->tend,0); CHKERRQ(ierr);
  ierr = PetscOptionsGetInt (0,"-timax",&w->timax,0); CHKERRQ(ierr);
  ierr = PetscOptionsGetInt (0,"-writeInterval",&w->writeInterval,0); CHKERRQ(ierr);

  ierr = PetscPrintf(comm, "DWorld Options\n"); CHKERRQ(ierr);
  ierr = PetscPrintf(comm, "CFL     = %e\n", w->CFL); CHKERRQ(ierr);
  ierr = PetscPrintf(comm, "dtmax   = %e\n", w->dtmax); CHKERRQ(ierr);
  ierr = PetscPrintf(comm, "dtframe = %e\n", w->dtframe); CHKERRQ(ierr);
  ierr = PetscPrintf(comm, "tend    = %e\n", w->tend); CHKERRQ(ierr);
  ierr = PetscPrintf(comm, "timax   = %d\n", w->timax); CHKERRQ(ierr);
  ierr = PetscPrintf(comm, "writeInterval   = %d\n", w->writeInterval); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
