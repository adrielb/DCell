#include "ImmersedInterfaceMethod.h"

struct _LocalCoor
{
  int Np;    // Max number of local coordinates
  int len;   // Current number of coordinates
  double *n; // Normal component
  double *s; // Tangential component
  double *r; // Tangential component
};

#define DROT drot_
void DROT(int *len, double* x, int *incx, double* y, int *incy, double *c, double *s); //from netlib blas
void cblas_drot(int len, double* x, int incx, double* y, int incy, double c, double s); //from intel mkl
double LocalCoor2DArcLengthFunction( double x, void *params);

#undef __FUNCT__
#define __FUNCT__ "LocalCoorCreate"
PetscErrorCode LocalCoorCreate( int Np, LocalCoor *lc )
{
  PetscErrorCode ierr;
  LocalCoor l;
  
  PetscFunctionBegin;

  if( Np < 5 ) 
    SETERRQ1( PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE,
        "Max points Np = %d need to be > 5", Np );
  
  ierr = PetscNew( struct _LocalCoor, &l); CHKERRQ(ierr);
  l->len = Np;
  l->Np = Np;
  ierr = PetscMalloc(Np*sizeof(double), &l->n); CHKERRQ(ierr);
  ierr = PetscMalloc(Np*sizeof(double), &l->s); CHKERRQ(ierr); 
  ierr = PetscMalloc(Np*sizeof(double), &l->r); CHKERRQ(ierr);
  
  *lc = l;
  
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "LocalCoorDestroy"
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
/*
void LocalCoor2DNormal( LocalCoor lc, LevelSet ls, IrregularNode *n )
{
  PetscReal nx, ny, h;
  
  nx = Bilinear2D(GridFunction2D_DerivX, ls->phi, ls->phi->d, n->x+n->ox, n->y+n->oy);
  ny = Bilinear2D(GridFunction2D_DerivY, ls->phi, ls->phi->d, n->x+n->ox, n->y+n->oy);
  h = sqrt( nx*nx + ny*ny );
  
  n->nx = nx / h;
  n->ny = ny / h;
  n->sx = ny /h;
  n->sy = -nx / h;
}
*/
void LocalCoor2DGetVecs( LocalCoor lc, PetscReal **s, PetscReal **n )
{
  *n = lc->n;
  *s = lc->s;
}

void LocalCoorSetLength( LocalCoor lc, PetscInt len)
{
  lc->len = len;
}

void LocalCoor2DSolve( LocalCoor lc, Coor dh, IrregularNode *n)
{
  int i, incx=1, incy=1;
  
  for( i = 0; i < lc->len; i++)
  {
    lc->s[i] = (lc->s[i] - n->X.x) * dh.x;
    lc->n[i] = (lc->n[i] - n->X.y) * dh.y;
  }
  
//  cblas_drot(lc->len, lc->s, 1, lc->n, 1, n->nx, n->ny);
  DROT(&lc->len, lc->s, &incx, lc->n, &incy, &n->nx, &n->ny);

  for( i = 0; i < lc->len; i++)
  {
    lc->s[i] = PetscSign(lc->s[i]) * PetscSqrtReal(
                ( lc->s[i] * lc->s[i] ) +
                ( lc->n[i] * lc->n[i] ) );
  }
}

void LocalCoor3DTangential( IrregularNode *n )
{
  PetscReal h;

  if( n->nx * n->nx + n->ny * n->ny >=
      n->nx * n->nx + n->nz * n->nz )
  {
    h = sqrt( n->nx * n->nx + n->ny * n->ny );
    n->sx =  n->ny / h;
    n->sy = -n->nx / h;
    n->sz = 0;

    n->rx =  n->nx * n->nz;
    n->ry =  n->ny * n->nz;
    n->rz = -n->nx * n->nx - n->ny * n->ny;
  } else {
    h = sqrt( n->nx * n->nx + n->nz * n->nz );
    n->sx = n->nz / h;
    n->sy = 0;
    n->sz =-n->nx / h;

    n->rx = -n->nx * n->ny;
    n->ry =  n->nx * n->nx + n->nz * n->nz;
    n->rz = -n->ny * n->nz;
  }
  h = sqrt( n->rx * n->rx + n->ry * n->ry + n->rz * n->rz );
  n->rx /= h;
  n->ry /= h;
  n->rz /= h;
}

void LocalCoor3DSolve( LocalCoor lc, Coor dh, IrregularNode *N )
{
  PetscReal *n = lc->n, *s = lc->s, *r = lc->r;
  PetscReal x, y, z;

  int i;
  for ( i = 0; i < lc->len; ++i)
  {
    x = (n[i] - N->X.x) * dh.x;
    y = (s[i] - N->X.y) * dh.y;
    z = (r[i] - N->X.z) * dh.y;
    n[i] = N->nx * x + N->ny * y + N->nz * z;
    s[i] = N->sx * x + N->sy * y + N->sz * z;
    r[i] = N->rx * x + N->ry * y + N->rz * z;
  }

  for ( i = 0; i < lc->len; ++i)
  {
    s[i] = PetscSign(s[i]) * PetscSqrtReal( s[i]*s[i] + n[i]*n[i] );
    r[i] = PetscSign(r[i]) * PetscSqrtReal( r[i]*r[i] + n[i]*n[i] );
  }
}

void LocalCoor3DGetVecs( LocalCoor lc, PetscReal **s, PetscReal **r, PetscReal **n )
{
  *n = lc->n;
  *s = lc->s;
  *r = lc->r;
}
