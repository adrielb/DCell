#include "ImmersedInterfaceMethod.h"
//#include "gsl/gsl_integration.h"
//#include "gsl/gsl_errno.h"

struct _LocalCoor
{
  int Np;    // Max number of local coordinates
  int len;   // Current number of coordinates
  double *n; // Normal component
  double *s; // Tangential component
  double *r; // Tangential component
//  double C,D;  // Hermite interpolation parameters
//  double epsabs, epsrel; // Integration tolerances
//  gsl_function F;
};


void drot_(int len, double* x, int incx, double* y, int incy, double c, double s);
void cblas_drot(int len, double* x, int incx, double* y, int incy, double c, double s);
double LocalCoor2DArcLengthFunction( double x, void *params);

#undef __FUNCT__
#define __FUNCT__ "LocalCoor2DCreate"
PetscErrorCode LocalCoorCreate( int Np, LocalCoor *lc )
{
  PetscErrorCode ierr;
  LocalCoor l;
  
  PetscFunctionBegin;
  
  if( Np < 5 ) 
    SETERRQ1( PETSC_ERR_ARG_OUTOFRANGE, 
        "Max points Np = %d need to be > 5", Np );
  
  ierr = PetscNew( struct _LocalCoor, &l); CHKERRQ(ierr);
  l->len = Np;
  l->Np = Np;
  ierr = PetscMalloc(Np*sizeof(double), &l->n); CHKERRQ(ierr);
  ierr = PetscMalloc(Np*sizeof(double), &l->s); CHKERRQ(ierr); 
  ierr = PetscMalloc(Np*sizeof(double), &l->r); CHKERRQ(ierr);
  
  /* gsl integration for arc length calculation
//  gsl_set_error_handler_off(); // prevent gsl from calling abort on error 
  l->epsabs = 0;
  l->epsrel = 1e-7;
//  l->F.function = &LocalCoor2DArcLengthFunction;
//  l->F.params = l;
  */
  
  *lc = l;
  
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "LocalCoor2DDestroy"
PetscErrorCode LocalCoorDestroy( LocalCoor lc )
{
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  ierr = PetscFree( lc->n ); CHKERRQ(ierr);
  ierr = PetscFree( lc->s ); CHKERRQ(ierr);
  ierr = PetscFree( lc->r ); CHKERRQ(ierr);
  ierr = PetscFree( lc    ); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

void LocalCoor2DNormal( LocalCoor lc, LevelSet ls, IrregularNode *n )
{
  PetscReal nx, ny, h;
  
  nx = Bilinear2D(GridFunction2D_DerivX, ls->g, n->x+n->ox, n->y+n->oy);
  ny = Bilinear2D(GridFunction2D_DerivY, ls->g, n->x+n->ox, n->y+n->oy);
  h = sqrt( nx*nx + ny*ny );
  
  n->nx = nx / h;
  n->ny = ny / h;
}

void LocalCoor2DGetVecs( LocalCoor lc, PetscReal **s, PetscReal **n )
{
  *n = lc->n;
  *s  = lc->s;
}

void LocalCoorSetLength( LocalCoor lc, PetscInt len)
{
  lc->len = len;
}

void LocalCoor2DSolve( LocalCoor lc, IrregularNode *n)
{
  PetscReal X = n->x+n->ox, Y = n->y+n->oy;
  int i;
  
  for( i = 0; i < lc->len; i++)
  {
    lc->s[i] -= X;
    lc->n[i] -= Y;
  }
  
//  cblas_drot(lc->len, lc->s, 1, lc->n, 1, n->nx, n->ny);
  drot_(lc->len, lc->s, 1, lc->n, 1, n->nx, n->ny);
}

/*
void LocalCoor2DSolveStencil( LocalCoor lc, IrregularNode *n )
{
  int i;
  for( i = 0; i < 5; ++i)
  {
    lc->s[i] = Stencil2Dx[i] - n->ox;
    lc->n[i]  = Stencil2Dy[i] - n->oy;
  }
  drot_(5, lc->s, 1, lc->n, 1, n->nx, n->ny);
}
*/

/*
int gsl_integration_qng (
    const gsl_function * f, 
    double a, double b, 
    double epsabs, double epsrel, 
    double * result, double * abserr, 
    size_t * neval);

void LocalCoor2DToArcLength( LocalCoor2D lc, IrregularNode *n,
    PetscReal n1, PetscReal n2, int i, PetscReal *s )
{
  PetscReal e = lc->eta[i], x = lc->xi[i];
  int err;
  size_t neval;
  double abserr;
  
  if( x == 0. )
  {
    *s = 0;
    return;
  }
  
  n1 = n->nx * n1 + n->ny * n2;
  n2 = n->nx * n2 - n->ny * n1;
  
  lc->C = -(x*n2 + 3*n1*e) / (x*x   * n1);
  lc->D =  (x*n2 + 2*n1*e) / (x*x*x * n1);
  
  err = gsl_integration_qng(
      &lc->F, 
      0, x, 
      lc->epsabs, lc->epsrel, 
      s, &abserr, 
      &neval);
  
  if( err ) // handel error 
  {
    //TODO: use petsc seterrq
    if( err == GSL_ETOL )
    {
      printf("err: %10.10f\t(%d,%d)\n", abserr, n->x, n->y);
      *s = PetscSign(x) * sqrt( e*e + x*x );
    }
    else
      exit(1);
  }

}
*/

/*
double LocalCoor2DArcLengthFunction( double x, void *params)
{
  LocalCoor lc = (LocalCoor) params;
  return sqrt( 1 + PetscSqr(2*lc->C*x + 3*lc->D*x*x) );
}
*/
