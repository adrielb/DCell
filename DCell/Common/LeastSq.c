//#include "mkl.h"
//#include "mkl_cblas.h"
#include "Common.h"
#define DGELSS dgelss_

void PrintMat( int numrows, int numcols, double *mat);
void DGELSS(int *m,int *n,int *nrhs,double *a,int *lda,double *b,int *ldb,
           double *s,double *rcond,int *rank,double *work,int *lwork,int *info);

struct _LeastSq {
  MPI_Comm comm;
  int Np;       // Maximum # of points
  int m;        // Current number of points
  int n;        // Number of basis (1,x,x^2) or (1,x,y,x^2,y^2,xy)
  double *a;    // [[...1...][...x...][...x^2...]]
  int nrhs;     // Input:  (m x 1) dependent variable to interpolate
  double *b;    // Output: (n x 1) [g, dg, ddg] 
  double rcond; // -1: all s.v. 
  int lwork;    // workspace size calculated during construction
  double *work; // workspace
  int rank;
  double *sv;   // Singular Values
  int info;
  double *s, *r;// Input: (m x 1) surface parameterization
};

#undef __FUNCT__
#define __FUNCT__ "LeastSqCreate"
PetscErrorCode LeastSqCreate( int Np, PetscBool is2D, LeastSq *ls )
{
  PetscErrorCode ierr;
  LeastSq lsq;
  
  PetscFunctionBegin;
  ierr = PetscNew(struct _LeastSq, &lsq); CHKERRQ(ierr);
  
  lsq->Np = Np;  
  lsq->m = Np;
  lsq->n = is2D ? 3 : 6;
  ierr = PetscMalloc(lsq->m*lsq->n*sizeof(double), &lsq->a); CHKERRQ(ierr);
  
  lsq->nrhs = 1;
  ierr = PetscMalloc( lsq->Np*lsq->nrhs*sizeof(double), &lsq->b); CHKERRQ(ierr);
  ierr = PetscMalloc( lsq->Np*lsq->nrhs*sizeof(double), &lsq->s); CHKERRQ(ierr);
  ierr = PetscMalloc( lsq->Np*lsq->nrhs*sizeof(double), &lsq->r); CHKERRQ(ierr);
  
  ierr = PetscMalloc( lsq->m*sizeof(double), &lsq->sv); CHKERRQ(ierr);
  lsq->rcond = -1;
  
  //Run dgelss to find lwork value
  PetscReal work;
  lsq->lwork = -1;
  DGELSS( &lsq->m, &lsq->n, &lsq->nrhs, lsq->a, &lsq->m,
    lsq->b, &lsq->m, lsq->sv, &lsq->rcond, &lsq->rank,
    &work, &lsq->lwork, &lsq->info ); CHKERRQ(lsq->info);
  lsq->lwork = work;
  
  ierr = PetscMalloc( lsq->lwork*sizeof(double), &lsq->work); CHKERRQ(ierr);
  
  *ls = lsq;
  
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "LeastSqDestroy"
PetscErrorCode LeastSqDestroy( LeastSq ls )
{
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  ierr = PetscFree( ls->a ); CHKERRQ(ierr);
  ierr = PetscFree( ls->b ); CHKERRQ(ierr);
  ierr = PetscFree( ls->s ); CHKERRQ(ierr);
  ierr = PetscFree( ls->r ); CHKERRQ(ierr);
  ierr = PetscFree( ls->sv ); CHKERRQ(ierr);
  ierr = PetscFree( ls->work ); CHKERRQ(ierr);
  ierr = PetscFree( ls ); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "LeastSqGetVecs"
PetscErrorCode LeastSqGetVecs( LeastSq ls, double **s, double **r, double **g, int *len )
{
  PetscFunctionBegin;
  if( s )   *s   = ls->s;
  if( r )   *r   = ls->r;
  if( g )   *g   = ls->b;
  if( len ) *len = ls->m;
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "LeastSqSetNumPoints"
PetscErrorCode LeastSqSetNumPoints( LeastSq ls, int m )
{
  PetscFunctionBegin;
  if( m > ls->Np ) 
    SETERRQ2( ls->comm, PETSC_ERR_ARG_OUTOFRANGE,
      "Max points: Np = %d < m = %d given.", ls->Np, m );
  if( m < ls->n ) 
    SETERRQ2( ls->comm, PETSC_ERR_ARG_OUTOFRANGE,
      "Min points: need n = %d > m = %d given.", ls->n, m );
  
  ls->m = m;
  
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "LeastSqSolve"
PetscErrorCode LeastSqSolve( LeastSq ls )
{
  double *s = ls->s, *r = ls->r;
  int m = ls->m;
  int i;
  
  PetscFunctionBegin;
  
//form the matrix 'a' from 's'
  if( ls->n == 3 ) {
    for( i = 0; i < ls->m; i++ )
    {
      ls->a[i + 0*m] = 1.0;
      ls->a[i + 1*m] = s[i];
      ls->a[i + 2*m] = 0.5 * s[i] * s[i];
    //    ls->a[0 + i*3] = one;
    //    ls->a[1 + i*3] = s[i];
    //    ls->a[2 + i*3] = 0.5 * s[i] * s[i];
    }
  } else {
    for( i = 0; i < ls->m; i++ )
    {
      ls->a[i + 0*m] = 1.0;
      ls->a[i + 1*m] = s[i];
      ls->a[i + 2*m] = r[i];
      ls->a[i + 3*m] = 0.5 * s[i] * s[i];
      ls->a[i + 4*m] = 0.5 * r[i] * r[i];
      ls->a[i + 5*m] = s[i] * r[i];
    }
  }
//  printf("\na\n");
//  PrintMat( ls->m, ls->Np, ls->a);
//  
//  printf("\nb\n");
//  PrintMat( ls->n, ls->m, ls->b);
    
  // [g dg ddg] = min | A x - y |
  DGELSS(
    &ls->m    /* num rows in A */,   // 7 
    &ls->n    /* num cols in A */,   // 3
    &ls->nrhs /* num cols in B */,   // 1
    ls->a     /* matrix A for SVD */, 
    &ls->m    /* first dim of A = max( 1, m )*/,  // 3 
    ls->b     /* matrix B of RHS (become solution matrix X) */, 
    &ls->m    /* first dim of B = max( 1, m, n ) */, // 7
    ls->sv    /* matrix of s.v */, 
    &ls->rcond/* relative condition s.v. treated as zero */, //-1 
    &ls->rank /*  */, 
    ls->work  /*  */, 
    &ls->lwork /* size of work */, 
    &ls->info /*  */
  );
  if( ls->info )
  {
    printf("1\tm:%d\n", ls->m);
    printf("2\tn:%d\n", ls->n);
    printf("3\tnrhs:%d\n", ls->nrhs);
    PrintMat( ls->m, ls->n, ls->a );
    SETERRQ(ls->comm, 1,"MKL ERROR in LeastSq");
  }
    
  PetscFunctionReturn(0);
}

void PrintMat( int numrows, int numcols, double *mat)
{
  int i, j;
  
  for( j = 0; j < numcols; j++ )
  {
    for( i = 0; i < numrows; i++ )
    {
//      printf("%d\n",i + numrows * j);
      printf( "%2.5f\t\t", mat[ i + numrows * j ] );
    }
    printf("\n");
  }
}


/*
  for( int i = 0; i < Np*nrhs; i++)
    printf("b %d: %f\n", i, b[i]);
  for( int i = 0; i < 3; i++)
    printf("s %d: %f\n", i, s[i]);
  for( int i = 0; i < m*n; i++)
    printf("a %d: %f\n", i, a[i]);
*/
