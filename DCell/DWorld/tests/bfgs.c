#include "DWorld.h"
#include "petscblaslapack.h"

void PrintVec( const char *name, PetscReal *v, int len);
void PrintMat( const char *name, PetscReal *m, int len);

// A := alpha*x*y' + A,
int dger_(int *m, int *n, double *alpha,
  double *x, int *incx, double *y, int *incy,
  double *a, int *lda);

void f( PetscReal *x, PetscReal *g ) {
  g[0] =  exp(x[0]-1) + 2*(x[0] - x[1]);
  g[1] = -exp(1-x[1]) - 2*(x[0] - x[1]);
}

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc, char **args) {

  PetscErrorCode ierr;
  PetscFunctionBegin;
  ierr = DCellInit(); CHKERRQ(ierr);


  int i, j, k;
  PetscReal zero = 0, one = 1, negone = -1;
  PetscReal sTy, yTBy, a1;
  PetscBLASInt n = 2;
  int lda = n;
  PetscReal *m1, *m2, *ssT, *BysT, *syTB, *B; // matrices
  PetscReal *s, *y, *d, *x0, *x1, *g0, *g1, *By; // vectors
  int inc = 1;
  const char nop = 'N';
  int maxIter = 8;
  PetscReal tol = 1e-6;
  PetscReal lambda = 1;

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

  // set [B] = [I]
  for (j = 0; j < n; ++j) {
    for (i = 0; i < n; ++i) {
      B[i+n*j] = i==j ? 1 : 0;
    }
  }
  PrintMat("B",B,n);

  // set phi0 = phi
  x0[0] = 0;
  x0[1] = 0;
  PrintVec( "x0", x0, n);

  for (i = 0; i < maxIter; ++i) {
    // Update RHS
    // Evaluate g0 = g( X = x0 )
    f( x0, g0 );
    PrintVec( "g0", g0, n);

    // Solve B.d = -g0
    // Mult  d = -B.g0
    BLASgemv_(&nop,&n,&n,&negone,B,&lda,g0,&inc,&zero,d,&inc);
    PrintVec( "d", d, n);

    // Update position
    // x1 = x0 + l d
    for (j = 0; j < n; ++j) {
      x1[j] = x0[j] + lambda*d[j];
    }
    PrintVec("x1", x1, n);

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
    printf("norm: %f\n", norm);

    // Evaluate g1 = g( X = x1 )
    f( x1, g1 );
    PrintVec("g1", g1,n);

    // Update Hessian
    // s = x1 - x0
    // y = g(x1) - g(x0)
    // B = B + (sT.y + yT.B.y)(s.sT) / (sT.y)^2 - ( B.y.sT + s.yT.B ) / (sT.y)
    for (j = 0; j < n; ++j) {
      s[j] = lambda*d[j];
      y[j] = g1[j] - g0[j];
    }
    PrintVec("s",s,n);
    PrintVec("y",y,n);

    // [ (sT.y + yT.B.y) / (sT.y)^2 ] (s.sT)
    BLASgemv_(&nop,&n,&n,&one,B,&lda,y,&inc,&zero,By,&inc); // By = B.y
    PrintVec("By",By,n);

    yTBy = BLASdot_(&n,y,&inc,By,&inc);  // yT.B.y = y.By
    sTy  = BLASdot_(&n,s,&inc,y,&inc); // sTy = sT.y
    a1 = ( sTy + yTBy ) / (sTy*sTy);
    printf("yTBy: %f\n", yTBy);
    printf("sTy: %f\n", sTy);
    printf("a1: %f\n", a1);

    PetscMemzero(ssT,matsize);
    PetscMemzero(m1,matsize);
    PetscMemzero(m2,matsize);
    PetscMemzero(BysT,matsize);
    PetscMemzero(syTB,matsize);
    dger_(&n,&n,&a1,s,&inc,s,&inc,ssT,&lda); // ssT = a1 * s.sT
    PrintMat("ssT",ssT,n);
    // ( B.y.sT + s.yT.B ) / (sT.y)
    dger_(&n,&n,&one,y,&inc,s,&inc,m1,&lda); // m1 = y.sT
    PrintMat("m1",m1,n);
    dger_(&n,&n,&one,s,&inc,y,&inc,m2,&lda); // m2 = s.yT
    PrintMat("m2",m2,n);
    dgemm_(&nop,&nop,&n,&n,&n,&one,B,&lda,m1,&lda,&zero,BysT,&lda); // BysT = B.m1
    PrintMat("BysT",BysT,n);
    dgemm_(&nop,&nop,&n,&n,&n,&one,m2,&lda,B,&lda,&zero,syTB,&lda); // syTB = m2.B
    PrintMat("syTB",syTB,n);

    for (j = 0; j < n*n; ++j) {
      B[j] = B[j] + ssT[j] - (BysT[j] + syTB[j])/sTy;
    }
    PrintMat("B",B,n);
    printf("=================================================\n");
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

  ierr = DCellFinalize(); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

void PrintVec( const char *name, PetscReal *v, int len)
{
  int i;
  printf("%s: ", name);
  for( i = 0; i < len; i++ )
    printf("%f\t",v[i]);
  printf("\n");
}

void PrintMat( const char *name, PetscReal *m, int len)
{
  int i,j;
  for( j = 0; j < len; j++ ) {
    printf("%s:\t", name);
    for( i = 0; i < len; i++ ) {
      printf("%f\t",m[j+len*i]);
    }
    printf("\n");
  }
  printf("\n");
}
