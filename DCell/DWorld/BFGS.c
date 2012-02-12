#include "DWorld.h"
#include "petscblaslapack.h"

// A := alpha*x*y' + A,
int dger_(int *m, int *n, double *alpha,
  double *x, int *incx, double *y, int *incy,
  double *a, int *lda);

PetscErrorCode DWorld_TimeStep( DWorld w );
PetscErrorCode DWorld_DebugWrite( DWorld w, int i );

#undef __FUNCT__
#define __FUNCT__ "DWorldSimulate_BFGS"
PetscErrorCode DWorldSimulate_BFGS(DWorld w)
{
  int i, j, c;
  PetscReal norm, g0norm, g1norm;
  PetscReal zero = 0, one = 1, negone = -1;
  PetscReal sTy, yTBy, a1;
  PetscBLASInt n;
  PetscReal *ssT, *BysT, *syTB, *B; // matrices
  PetscReal *s, *y, *d, *g0, *g1, *temp, *By, *yB; // vectors
  int inc = 1;
  const char nop = 'N', T = 'T';
  int lda;
  int maxIter = 30;
  PetscReal tolx = 0.0;
  PetscReal tolg = 0.1;
  PetscReal lambda = 1;
  FluidField fluid = w->fluid;
  size_t vecsize, matsize;
  DCellsArray dcells = w->dcells;
  PetscErrorCode ierr;

  PetscFunctionBegin;
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

  // Initialize start and Hessian
  // set [B] = [I]
  for (j = 0; j < n; ++j) {
    for (i = 0; i < n; ++i) {
      B[i+n*j] = i==j ? 1 : 0;
    }
  }

  // Update RHS
  // Evaluate g0 = g( X = x0 )
  // g(X) = X - Xn - dt*U(X)
  ierr = VecZeroEntries(fluid->rhs); CHKERRQ(ierr);
  ierr = DCellsArrayUpdateFluidFieldRHS( dcells, w->iim, fluid, w->t ); CHKERRQ(ierr);
  ierr = FluidFieldSolve( fluid ); CHKERRQ(ierr);
  ierr = DWorld_TimeStep( w ); CHKERRQ(ierr);
  ierr = DCellsArrayAdvectImplicitRHS( dcells, w->fluid, w->dt, g0 ); CHKERRQ(ierr);
  ierr = DWorld_DebugWrite( w, 0); CHKERRQ(ierr);

  // Check root
  //  | g(x0) |  < tolg
  g0norm = 0;
  for (j = 0; j < n; ++j)
    g0norm += g0[j]*g0[j];
  g0norm = sqrt(g0norm);
  ierr = PetscInfo1(0,"BFGS root norm: %f\n", g0norm); CHKERRQ(ierr);

  for (i = 0; i < maxIter; ++i) {
    // Solve B.d = -g0
    // Mult  d = -B.g0
    BLASgemv_(&nop,&n,&n,&negone,B,&lda,g0,&inc,&zero,d,&inc);

    lambda = 1;
    for (c = 0; c < 30; ++c) {
      // Update position
      // x1 = x0 + l d
//      ierr = DWorld_DebugWrite( w, 1); CHKERRQ(ierr);
      ierr = DCellsArrayAdvectImplicitUpdate( dcells, lambda, d); CHKERRQ(ierr);
//      ierr = DWorld_DebugWrite( w, i+1); CHKERRQ(ierr);
//      ierr = DWorld_DebugWrite( w, 2); CHKERRQ(ierr);

      int m;
      for (m = 0; m < 10; ++m) {
        // Evaluate g1 = g( X = x1 )
        // g(X) = X - Xn - dt*U(X)
        ierr = VecZeroEntries(fluid->rhs); CHKERRQ(ierr);
        ierr = DCellsArrayUpdateFluidFieldRHS( dcells, w->iim, fluid, w->t ); CHKERRQ(ierr);
        ierr = FluidFieldSolve( fluid ); CHKERRQ(ierr);
        ierr = DCellsArrayAdvectImplicitRHS( dcells, w->fluid, w->dt, g1 ); CHKERRQ(ierr);

        // Check root
        //  | g(x1) |  < tolg
        g1norm = 0;
        for (j = 0; j < n; ++j) {
          g1norm += g1[j]*g1[j];
        }
        g1norm = sqrt(g1norm);
        if( g1norm < tolg )
          break;
      }

      if( g1norm > g0norm ) {
        // reset psi back to x0
        ierr = DCellsArrayAdvectImplicitUpdate( dcells, -lambda, d); CHKERRQ(ierr);
//        ierr = DWorld_DebugWrite( w, 3); CHKERRQ(ierr);
//        lambda = lambda / 2.;
      } else {
        g0norm = g1norm;
        break;
      }
    }
    ierr = PetscInfo1(0,"BFGS root norm: %f\n", g1norm); CHKERRQ(ierr);

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
      B[j] = B[j] + ssT[j] - (BysT[j] + syTB[j])/sTy;
    }
    PetscPrintf(PETSC_COMM_WORLD,"BFGS iter: %d\n",i);

    // Update g0 <-- g1
    temp = g0;
    g0 = g1;
    g1 = temp;
  } // i < maxIter

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
  PetscFunctionReturn(0);
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
  ierr = GridWrite(ls->psi->phi,i); CHKERRQ(ierr);
  ierr = FluidFieldWrite(w->fluid, i); CHKERRQ(ierr);

  ierr = LevelSetWriteIrregularNodeList(ls->psi, i); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
