#include "DWorld.h"
#include "petscblaslapack.h"

// max_i ( x )
int idamax_(int *n, double *dx, int *incx);

PetscErrorCode DWorld_TimeStep(DWorld w);
PetscErrorCode DWorld_DebugWrite(DWorld w, int i);
PetscErrorCode DWorldSimulate_ImplicitStep(DWorld w);
PetscReal gnorm(int n, PetscReal *g);

#undef __FUNCT__
#define __FUNCT__ "DWorldSimulate_Implicit"
PetscErrorCode DWorldSimulate_Implicit(DWorld w) {
  PetscErrorCode ierr;

  PetscFunctionBegin;
  if( w->g0array == PETSC_NULL ) {
    ierr = ArrayCreate("g0array",sizeof(double),16,&w->g0array); CHKERRQ(ierr);
    ierr = ArrayCreate("g1array",sizeof(double),16,&w->g1array); CHKERRQ(ierr);
    ierr = ArrayCreate("darray", sizeof(double),16,&w->darray); CHKERRQ(ierr);
  }

  // defaults
  w->MAX_STEPS_PICARD = 15;
  w->MAX_STEPS_LINE = 3;
  w->MAX_STEPS_DT = 10;
  w->MAX_D = 0.1;
  w->TOL_G = 0.01;

  // cmd line opts
  ierr = PetscOptionsGetInt( 0, "-max_steps_picard", &w->MAX_STEPS_PICARD, 0); CHKERRQ(ierr);
  ierr = PetscOptionsGetInt( 0, "-max_steps_line", &w->MAX_STEPS_LINE, 0); CHKERRQ(ierr);
  ierr = PetscOptionsGetInt( 0, "-max_steps_dt", &w->MAX_STEPS_DT, 0); CHKERRQ(ierr);
  ierr = PetscOptionsGetReal( 0, "-implicit_max_d", &w->MAX_D, 0); CHKERRQ(ierr);
  ierr = PetscOptionsGetReal( 0, "-implicit_tol_g", &w->TOL_G, 0); CHKERRQ(ierr);

  for (w->ti = 0; w->t < w->tend && w->ti < w->timax; ++w->ti) {
    ierr = DWorldSimulate_ImplicitStep(w); CHKERRQ(ierr);
    ierr = DWorldWrite(w, w->ti); CHKERRQ(ierr);
    ierr = DWorldPrintStep(w); CHKERRQ(ierr);
    w->t = w->t + w->dt;
  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DWorldSimulate_ImplicitStep"
PetscErrorCode DWorldSimulate_ImplicitStep(DWorld w) {
  int p, c, s, j;
  PetscReal g0norm, g1norm = PETSC_MAX_REAL;
  PetscBLASInt n;
  PetscReal *d, *g0, *g1;
  PetscReal *temp; // swaps g0 <-> g1
  int inc = 1;
  PetscReal lambda;
  FluidField fluid = w->fluid;
  DCellsArray dcells = w->dcells;
  PetscLogDouble t1, t2;
  int count = 0;
  const int MAX_STEPS_PICARD = w->MAX_STEPS_PICARD;
  const int MAX_STEPS_LINE = w->MAX_STEPS_LINE;
  const int MAX_STEPS_DT = w->MAX_STEPS_DT;
  const PetscReal TOL_G = w->TOL_G;
  const PetscReal MAX_D = w->MAX_D;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscGetTime(&t1); CHKERRQ(ierr);
  // set psi_0 = phi_n
  n = 0;
  ierr = DCellsArrayAdvectImplicitInit(dcells, &n); CHKERRQ(ierr);
  ierr = ArraySetSize( w->g0array, n); CHKERRQ(ierr);
  ierr = ArraySetSize( w->g1array, n); CHKERRQ(ierr);
  ierr = ArraySetSize( w->darray, n); CHKERRQ(ierr);
  g0   = ArrayGetData( w->g0array );
  g1   = ArrayGetData( w->g1array );
  d    = ArrayGetData( w->darray );

  // Update RHS
  // Evaluate g0 = g( X = x0 )
  // g(X) = X - Xn - dt*U(X)
  ierr = VecZeroEntries(fluid->rhs); CHKERRQ(ierr);
  ierr = DCellsArrayUpdateFluidFieldRHS(dcells, w->iim, fluid, w->t); CHKERRQ(ierr);
  ierr = FluidFieldSolve(fluid); CHKERRQ(ierr);
  ierr = VecCopy(fluid->vel, fluid->vel0); CHKERRQ(ierr);
  ierr = DWorld_TimeStep(w); CHKERRQ(ierr);
  ierr = DCellsArrayAdvectImplicitRHS(dcells, w->fluid, w->dt, g0); CHKERRQ(ierr);
  g0norm = gnorm(n, g0);
  ierr = PetscInfo1(0,"root norm = %f\n", g0norm); CHKERRQ(ierr);
  if (g0norm < TOL_G) {
    ierr = PetscInfo2(0,"Converged: root norm = %f < tolg = %f\n", g0norm, TOL_G); CHKERRQ(ierr);
    goto converged;
  }

  for (s = 0; s < MAX_STEPS_DT; ++s) {
    for (p = 0; p < MAX_STEPS_PICARD; ++p) {
      ierr = BernoulliDist( rnd, n, &delta); CHKERRQ(ierr);

      // g0 = L( q + c d )
      ierr = DCellsArrayAdvectImplicitUpdate(dcells, ck, delta[j] ); CHKERRQ(ierr);
      ierr = LossFunction( w, g0 ); CHKERRQ(ierr);
      ierr = DCellsArrayAdvectImplicitUpdate(dcells, -ck, delta[j] ); CHKERRQ(ierr);

      // g1 = L( q + c d - 2 c d == q - c d )
      ierr = DCellsArrayAdvectImplicitUpdate(dcells, -ck, delta[j] ); CHKERRQ(ierr);
      ierr = LossFunction( w, g1 ); CHKERRQ(ierr);
      ierr = DCellsArrayAdvectImplicitUpdate(dcells, ck, delta[j] ); CHKERRQ(ierr);

      for (j = 0; j < n; ++j) {
        d[j] = g1[j] - g0[j] / (2 * ck * delta[j] );
      }

      lambda = 1;
      for (c = 0; c < MAX_STEPS_LINE; ++c) {
        //* Check step size isn't too large
        int maxdi = idamax_(&n, d, &inc);
        PetscReal maxd = PetscAbs(d[maxdi]*lambda);
        if (maxd > MAX_D) {
          lambda = lambda * MAX_D / maxd;
        }
        ierr = PetscInfo1(0,"max lambda*d[j] = %f\n", lambda * maxd); CHKERRQ(ierr);

        // Update position
        // x1 = x0 + l d
        ierr = DCellsArrayAdvectImplicitUpdate(dcells, -lambda, d); CHKERRQ(ierr);



        g1norm = gnorm(n, g1);
        ierr = PetscInfo1(0,"tentative root norm = %f\n", g1norm); CHKERRQ(ierr);

        /*
        if( w->ti == 4 ) {
          ierr = DWorld_DebugWrite( w, count); CHKERRQ(ierr);
          count++;
        } else {
          if( w->ti == 5 )
            exit(1);
        }
        */

        if (g1norm < g0norm) {
          g0norm = g1norm;
          break;
        }

        reset:
        // reset psi back to x0
        ierr = DCellsArrayAdvectImplicitUpdate(dcells, lambda, d); CHKERRQ(ierr);
        lambda = lambda / 2.;
      } // c < MAX_STEP_ITER

      if (c == MAX_STEPS_LINE) {
        ierr = PetscInfo1(0,"MAX_STEPS_LINE = %d reached!\n", MAX_STEPS_LINE); CHKERRQ(ierr);
        break;
      }

      ierr = PetscInfo1(0,"root norm = %f\n", g1norm); CHKERRQ(ierr);

      // Check root
      //  | g(x1) |  < tolg
      if (g1norm < TOL_G) {
        ierr = PetscInfo2(0,"Converged: root norm = %f < tolg = %f\n", g1norm, TOL_G); CHKERRQ(ierr);
        goto converged;
      }

      // Update g0 <-- g1
      temp = g0;
      g0 = g1;
      g1 = temp;
    } // p < MAX_STEPS_PICARD

    if (p == MAX_STEPS_PICARD) {
      ierr = PetscInfo1( 0, "MAX_STEPS_PICARD = %d reached! \n", MAX_STEPS_PICARD); CHKERRQ(ierr);
    }

    ierr = DWorldSetTimeStep(w, w->dt / 2.); CHKERRQ(ierr);
    ierr = VecCopy(fluid->vel0, fluid->vel); CHKERRQ(ierr);
    ierr = DCellsArrayAdvectImplicitRHS(dcells, w->fluid, w->dt, g0); CHKERRQ(ierr);
    g0norm = gnorm(n, g0);
    ierr = PetscInfo1(0,"root norm = %f\n", g0norm); CHKERRQ(ierr);
    if (g0norm < TOL_G) {
      ierr = PetscInfo2(0,"Converged: root norm = %f < tolg = %f\n", g0norm, TOL_G); CHKERRQ(ierr);
      goto converged;
    }
  } // s < MAX_DT_STEPS

  if (s == MAX_STEPS_DT) {
    ierr = PetscInfo1( 0, "MAX_STEPS_DT = %d reached! \n", MAX_STEPS_DT); CHKERRQ(ierr);
    ierr = DCellFinalize(); CHKERRQ(ierr);
    exit(1);
  }

  converged:
  ierr = DCellsArrayAdvectImplicitReinit(dcells, w->dt); CHKERRQ(ierr);
  ierr = PetscGetTime(&t2); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "TS Solve: %f sec\n", t2 - t1); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "BernoulliDist"
PetscErrorCode BernoulliDist( PetscRandom rnd, int n, PetscReal *dist ) {
  int i;
  PetscReal val;
  PetscErrorCode ierr;
  PetscFunctionBegin;
  for (i = 0; i < n; ++i) {
    ierr = PetscRandomGetValueReal(rnd,&val); CHKERRQ(ierr);
    dist[i] = val < 0.5 ? -1 : 1;
  }
  PetscFunctionReturn(0);
}

PetscReal gnorm(int n, PetscReal *g) {
  //* max norm
  int j;
  PetscReal abs, norm = 0;
  for (j = 0; j < n; ++j) {
    abs = PetscAbs(g[j]);
    norm = norm < abs ? abs : norm;
  }
  return norm;
  //*/

  /* 2-norm
   int inc = 1;
   return dnrm2_( &n, g, &inc);
   */
}

#undef __FUNCT__
#define __FUNCT__ "LossFunction"
PetscErrorCode LossFunction(DWorld w, PetscReal *g) {
  FluidField fluid = w->fluid;
  DCellsArray dcells;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  // Evaluate g1 = g( X = x1 )
  // g(X) = X - Xn - dt*U(X)
  ierr = VecZeroEntries(fluid->rhs); CHKERRQ(ierr);
  ierr = DCellsArrayUpdateFluidFieldRHS(dcells, w->iim, fluid, w->t); CHKERRQ(ierr);
  ierr = FluidFieldSolve(fluid); CHKERRQ(ierr);
  ierr = VecAXPBY(fluid->vel, 0.5, 0.5, fluid->vel0); CHKERRQ(ierr);
  PetscReal maxVel, dtcfl;
  ierr = FluidFieldMaxVelocityMag(w->fluid, &maxVel); CHKERRQ(ierr);
  dtcfl = w->CFL * w->fluid->dh.x / maxVel;
  if (dtcfl < w->dt) {
    // velocity too fast for current cfl condition
    ierr = PetscInfo2(0,"Computed velocity too fast for current CFL condition: dtcfl = %f < dt = %f\n", dtcfl, w->dt); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

  ierr = DCellsArrayAdvectImplicitRHS(dcells, w->fluid, w->dt, g); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DWorld_TimeStep"
PetscErrorCode DWorld_TimeStep(DWorld w) {
  PetscReal dt;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = FluidFieldMaxVelocityMag(w->fluid, &w->maxVel); CHKERRQ(ierr);
  w->dtcfl = w->CFL * w->fluid->dh.x / w->maxVel;
  dt = PetscMin(w->dtmax,w->dtcfl);
  ierr = DWorldSetTimeStep(w, dt); CHKERRQ(ierr);
  printf("dt: %f\n", dt);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DWorld_DebugWrite"
PetscErrorCode DWorld_DebugWrite(DWorld w, int i) {
  DCell dcell;
  LevelSet ls;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = ArrayGetP(w->dcells->dcells, 0, &dcell); CHKERRQ(ierr);
  ls = dcell->lsPlasmaMembrane;
  ierr = GridWrite(ls->phi, i); CHKERRQ(ierr);
  ierr = GridWrite(ls->tmp, i); CHKERRQ(ierr);
  ierr = GridWrite(ls->psi->phi, i); CHKERRQ(ierr);
  ierr = FluidFieldWrite(w->fluid, i); CHKERRQ(ierr);

  ierr = LevelSetUpdateIrregularNodeList(ls); CHKERRQ(ierr);
  ierr = LevelSetWriteIrregularNodeList(ls, i); CHKERRQ(ierr);
  ierr = LevelSetWriteIrregularNodeList(ls->psi, i); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
