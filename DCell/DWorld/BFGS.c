#include "DWorld.h"
#include "petscblaslapack.h"

// A := alpha*x*y' + A,
int dger_(int *m, int *n, double *alpha,
  double *x, int *incx, double *y, int *incy,
  double *a, int *lda);

// max_i ( x )
int idamax_(int *n, double *dx, int *incx);

PetscErrorCode DWorld_TimeStep( DWorld w );
PetscErrorCode DWorld_DebugWrite( DWorld w, int i );
PetscErrorCode DWorldSimulate_BFGSstep(DWorld w);
void InitHessian( int n, PetscReal *B );
PetscReal gnorm( int n, PetscReal *g);

#undef __FUNCT__
#define __FUNCT__ "DWorldSimulate_BFGS"
PetscErrorCode DWorldSimulate_BFGS(DWorld w)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  for ( w->ti = 0; w->t < w->tend && w->ti < w->timax; ++w->ti) {
    ierr = DWorldSimulate_BFGSstep( w ); CHKERRQ(ierr);
    ierr = DWorldWrite( w, w->ti); CHKERRQ(ierr);
    ierr = DWorldPrintStep(w); CHKERRQ(ierr);
    w->t = w->t + w->dt;
  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DWorldSimulate_BFGSstep"
PetscErrorCode DWorldSimulate_BFGSstep(DWorld w)
{
  int i, j, c;
  PetscReal g0norm, g1norm;
  PetscReal zero = 0, one = 1, negone = -1;
  PetscReal sTy, yTBy, a1;
  PetscBLASInt n;
  PetscReal *ssT, *BysT, *syTB, *B; // matrices (4)
  PetscReal *s, *y, *d, *g0, *g1, *By, *yB; // vectors (7)
  PetscReal *temp; // swaps g0 <-> g1
  int inc = 1;
  const char nop = 'N', T = 'T';
  int lda;
  int MAX_BFGS_ITER = 15;
  int MAX_STEP_ITER = 3;
  PetscReal tolg = 0.01;
  PetscReal MAX_D = 0.05;
  PetscReal lambda = 1;
  FluidField fluid = w->fluid;
  size_t vecsize, matsize;
  DCellsArray dcells = w->dcells;
  PetscLogDouble t1, t2;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscGetTime(&t1); CHKERRQ(ierr);
  // set psi_0 = phi_n
  n = 0;
  ierr = DCellsArrayAdvectImplicitInit(dcells,&n); CHKERRQ(ierr);

  lda = n;
  vecsize = sizeof(double)*n;
  matsize = sizeof(double)*n*n;
  ierr = PetscMalloc(vecsize,&s); CHKERRQ(ierr);
  ierr = PetscMalloc(vecsize,&y); CHKERRQ(ierr);
  ierr = PetscMalloc(vecsize,&d); CHKERRQ(ierr);
  ierr = PetscMalloc(vecsize,&g0); CHKERRQ(ierr);
  ierr = PetscMalloc(vecsize,&g1); CHKERRQ(ierr);
  ierr = PetscMalloc(vecsize,&By); CHKERRQ(ierr);
  ierr = PetscMalloc(vecsize,&yB); CHKERRQ(ierr);
  ierr = PetscMalloc(matsize,&ssT); CHKERRQ(ierr);
  ierr = PetscMalloc(matsize,&BysT); CHKERRQ(ierr);
  ierr = PetscMalloc(matsize,&syTB); CHKERRQ(ierr);
  ierr = PetscMalloc(matsize,&B); CHKERRQ(ierr);

/*
  if( w->bfgs == NULL ) {
    ierr = ArrayCreate("bfgs", sizeof(double), n, &w->bfgs); CHKERRQ(ierr);
  }

  ierr = ArraySetSize(w->bfgs, 7*n + 4*n*n ); CHKERRQ(ierr);
  ierr = ArrayZero(w->bfgs); CHKERRQ(ierr);
  PetscReal *work = ArrayGetData(w->bfgs);
  s  = &work[0*n];
  y  = &work[1*n];
  d  = &work[2*n];
  g0 = &work[3*n];
  g1 = &work[4*n];
  By = &work[5*n];
  yB = &work[6*n];
  ssT = &work[7*n + 0*n*n];
  BysT= &work[7*n + 1*n*n];
  syTB= &work[7*n + 2*n*n];
  B   = &work[7*n + 3*n*n];
*/

  // Initialize Hessian
  // set [B] = [I]
  InitHessian(n,B);

  // Update RHS
  // Evaluate g0 = g( X = x0 )
  // g(X) = X - Xn - dt*U(X)
  ierr = VecZeroEntries(fluid->rhs); CHKERRQ(ierr);
  ierr = DCellsArrayUpdateFluidFieldRHS( dcells, w->iim, fluid, w->t ); CHKERRQ(ierr);
  ierr = FluidFieldSolve( fluid ); CHKERRQ(ierr);
  ierr = VecCopy(fluid->vel,fluid->vel0); CHKERRQ(ierr);
  ierr = DWorld_TimeStep( w ); CHKERRQ(ierr);
  ierr = DCellsArrayAdvectImplicitRHS( dcells, w->fluid, w->dt, g0 ); CHKERRQ(ierr);
//  ierr = DWorld_DebugWrite( w, 0); CHKERRQ(ierr);

  g0norm = gnorm(n, g0);
  ierr = PetscInfo1(0,"root norm = %f\n", g0norm); CHKERRQ(ierr);
int count = 0;
  for (i = 0; i < MAX_BFGS_ITER; ++i) {
    // Solve B.d = -g0
    // Mult  d = -B.g0
    BLASgemv_(&nop,&n,&n,&negone,B,&lda,g0,&inc,&zero,d,&inc);

    lambda = 1;
    for (c = 0; c < MAX_STEP_ITER; ++c) {
      //* Check step size isn't too large
      int maxdi = idamax_( &n, d, &inc);
      PetscReal maxd = PetscAbs(d[maxdi]*lambda);
      if( maxd > MAX_D ) {
        lambda = lambda * MAX_D / maxd;
      }
      ierr = PetscInfo1(0,"max lambda*d[j] = %f\n", lambda * maxd); CHKERRQ(ierr);

      // Update position
      // x1 = x0 + l d
      ierr = DCellsArrayAdvectImplicitUpdate( dcells, lambda, d); CHKERRQ(ierr);

      // Evaluate g1 = g( X = x1 )
      // g(X) = X - Xn - dt*U(X)
      ierr = VecZeroEntries(fluid->rhs); CHKERRQ(ierr);
      ierr = DCellsArrayUpdateFluidFieldRHS( dcells, w->iim, fluid, w->t ); CHKERRQ(ierr);
      ierr = FluidFieldSolve( fluid ); CHKERRQ(ierr);
      ierr = VecAXPBY(fluid->vel,0.5,0.5,fluid->vel0); CHKERRQ(ierr);
      ierr = DCellsArrayAdvectImplicitRHS( dcells, w->fluid, w->dt, g1 ); CHKERRQ(ierr);

      g1norm = gnorm( n, g1);
      ierr = PetscInfo1(0,"tentative root norm = %f\n", g1norm); CHKERRQ(ierr);

      if( w->ti == 10 ) {
        count++;
        ierr = DWorld_DebugWrite( w, count); CHKERRQ(ierr);
      }

      if( g1norm > g0norm ) {
        // reset psi back to x0
        ierr = DCellsArrayAdvectImplicitUpdate( dcells, -lambda, d); CHKERRQ(ierr);
        lambda = lambda / 2.;
      } else {
        g0norm = g1norm;
        break;
      }
    } // c < MAX_STEP_ITER

    if( c == MAX_STEP_ITER ) {
      ierr = PetscInfo1(0,"MAX_STEP_ITER = %d reached!\n", MAX_STEP_ITER); CHKERRQ(ierr);
      g0norm = g1norm;
//      exit(1);
//      continue;
    }

    ierr = PetscInfo1(0,"root norm = %f\n", g1norm); CHKERRQ(ierr);

    // Check root
    //  | g(x1) |  < tolg
    if( g1norm < tolg ) {
      ierr = PetscInfo2(0,"Converged: root norm = %f < tolg = %f\n", g1norm, tolg); CHKERRQ(ierr);
      break;
    }

    temp = g0;
    g0 = g1;
    g1 = temp;
    continue;

    // Update Hessian
    // s = x1 - x0 = l * d
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

    PetscMemzero(ssT,matsize);
    PetscMemzero(BysT,matsize);
    PetscMemzero(syTB,matsize);

    // ssT = a1 * s.sT
    dger_(&n,&n,&a1,s,&inc,s,&inc,ssT,&lda);
    // ( B.y.sT + s.yT.B ) / (sT.y)
    // yB = yT.B
    BLASgemv_(&T,&n,&n,&one,B,&lda,y,&inc,&zero,yB,&inc);
    dger_(&n,&n,&one,By,&inc,s,&inc,BysT,&lda); // BysT = By.sT
    dger_(&n,&n,&one,s,&inc,yB,&inc,syTB,&lda); // syTB = s.yB

    for (j = 0; j < n*n; ++j) {
      B[j] = B[j] + (ssT[j] - (BysT[j] + syTB[j])/sTy);
    }

    // Update g0 <-- g1
    temp = g0;
    g0 = g1;
    g1 = temp;
  } // i < MAX_BFGS_ITER

  if( i == MAX_BFGS_ITER ) {
    ierr = PetscInfo1( 0, "MAX_BFGS_ITER = %d reached! \n", MAX_BFGS_ITER); CHKERRQ(ierr);
    exit(1);
  }

  ierr = DCellsArrayAdvectImplicitReinit( dcells, w->dt); CHKERRQ(ierr);

  ierr = PetscFree(s); CHKERRQ(ierr);
  ierr = PetscFree(y); CHKERRQ(ierr);
  ierr = PetscFree(d); CHKERRQ(ierr);
  ierr = PetscFree(g0); CHKERRQ(ierr);
  ierr = PetscFree(g1); CHKERRQ(ierr);
  ierr = PetscFree(By); CHKERRQ(ierr);
  ierr = PetscFree(yB); CHKERRQ(ierr);
  ierr = PetscFree(ssT); CHKERRQ(ierr);
  ierr = PetscFree(BysT); CHKERRQ(ierr);
  ierr = PetscFree(syTB); CHKERRQ(ierr);
  ierr = PetscFree(B); CHKERRQ(ierr);

  ierr = PetscGetTime(&t2); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"BFGS Solve: %f sec\n", t2-t1); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

void InitHessian( int n, PetscReal *B )
{
  int i, j;
  for (j = 0; j < n; ++j) {
    for (i = 0; i < n; ++i) {
      B[i+n*j] = i==j ? 1 : 0;
    }
  }
}

PetscReal gnorm( int n, PetscReal *g) {
  //* max norm
  int j;
  PetscReal abs, norm;
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
#define __FUNCT__ "DWorld_TimeStep"
PetscErrorCode DWorld_TimeStep( DWorld w )
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = FluidFieldMaxVelocityMag( w->fluid, &w->maxVel); CHKERRQ(ierr);
  w->dtcfl = w->CFL * w->fluid->dh.x / w->maxVel;
  w->dt = PetscMin(w->dtmax,w->dtcfl);
  printf("dt: %f\n", w->dt);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DWorld_DebugWrite"
PetscErrorCode DWorld_DebugWrite( DWorld w, int i )
{
  DCell dcell;
  LevelSet ls;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = ArrayGetP(w->dcells->dcells,0,&dcell); CHKERRQ(ierr);
  ls = dcell->lsPlasmaMembrane;
  ierr = GridWrite(ls->phi,i); CHKERRQ(ierr);
  ierr = GridWrite(ls->tmp,i); CHKERRQ(ierr);
  ierr = GridWrite(ls->psi->phi,i); CHKERRQ(ierr);
  ierr = FluidFieldWrite(w->fluid, i); CHKERRQ(ierr);

  ierr = LevelSetUpdateIrregularNodeList(ls); CHKERRQ(ierr);
  ierr = LevelSetWriteIrregularNodeList(ls, i); CHKERRQ(ierr);
  ierr = LevelSetWriteIrregularNodeList(ls->psi, i); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


/* Check convergence
// | x1 - x0 | == | l*d | < tolx
norm = 0;
for (j = 0; j < n; ++j) {
  norm += lambda*d[j]*lambda*d[j];
}
norm = sqrt(norm);
if( norm < tolx ) {
  ierr = PetscInfo2(0,"BFGS:   dx norm = %f < tolx = %f converged.\n", norm, tolx ); CHKERRQ(ierr);
  break;
}
*/
