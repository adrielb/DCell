#include "PoissonSolver.h"
#include "Utilities.h"
#include "Grid.h"
#include "mkl_dfti.h"
#include "mkl_poisson.h"
#include <stdio.h>
#include <omp.h>

#undef __FUNCT__
#define __FUNCT__ "SpectralMethod2D"
PetscInt EVENT_SpectralMethod2D;
PetscErrorCode SpectralMethod2D(Grid g)
{
  PetscErrorCode ierr;
  MKL_INT nx=g->n.x-1, ny=g->n.y-1;
  /* Note that the size of the transform nx must be even !!! */
  MKL_INT ix, iy, i, stat;
  MKL_INT ipar[128];
  double ax, bx, ay, by;
  double *dpar, *f, *bd_ax, *bd_bx, *bd_ay, *bd_by;
  double q=0;
  PetscLogDouble p1, p2;
  DFTI_DESCRIPTOR_HANDLE xhandle = 0;
  char *BCtype;
    
  PetscFunctionBegin;
  PetscLogEventBegin(EVENT_SpectralMethod2D,0,0,0,0);

  dpar=(double*)malloc((13*nx/2+7)*sizeof(double));
  f=g->v1;
  bd_ax=(double*)malloc((ny+1)*sizeof(double));
  bd_bx=(double*)malloc((ny+1)*sizeof(double));
  bd_ay=(double*)malloc((nx+1)*sizeof(double));
  bd_by=(double*)malloc((nx+1)*sizeof(double));

  /* Defining the rectangular domain 0<x<1, 0<y<1 for 2D Laplace Solver */
  ax=0.0E0;  bx=1.0E0;
  ay=0.0E0;  by=1.0E0;

  BCtype = "NNNN";
  
  for(iy=0;iy<=ny;iy++)
  {
    bd_ax[iy]=0.;
    bd_bx[iy]=0.;
  }
  
  for(ix=0;ix<=nx;ix++)
  {
    bd_ay[ix]=0.;
    bd_by[ix]=0.;
  }

  /* Initializing ipar array to make it free from garbage */
  PetscMemzero(ipar, 128*sizeof(MKL_INT));

  d_init_Helmholtz_2D(&ax, &bx, &ay, &by, &nx, &ny, 
        BCtype, &q, ipar, dpar, &stat);CHKERRQ(stat);
  
  //Commit will alter RHS "f" and boundary values    
  d_commit_Helmholtz_2D(f, bd_ax, bd_bx, bd_ay, bd_by, 
        &xhandle, ipar, dpar, &stat); CHKERRQ(stat);

  d_Helmholtz_2D(f, bd_ax, bd_bx, bd_ay, bd_by, 
        &xhandle, ipar, dpar, &stat); CHKERRQ(stat);
  
  free_Helmholtz_2D(&xhandle, ipar, &stat); CHKERRQ(stat);
  
  free( dpar );
  free( bd_ax);
  free( bd_ay);
  free( bd_bx);
  free( bd_by);
  
  PetscLogEventEnd(EVENT_SpectralMethod2D,0,0,0,0);
  PetscFunctionReturn(0);
}

PetscErrorCode SpectralMethodSingle( int nn, PetscLogDouble *time)
{
  PetscErrorCode ierr;
  PetscLogDouble p1, p2;
  
  MKL_INT nx=nn, ny=nn, nz=nn;
  MKL_INT ix, iy, iz, i, stat;
  MKL_INT ipar[128];
  float ax, bx, ay, by, az, bz, lx, ly, lz, hx, hy, hz, xi, yi, zi, cx, cy, cz;
  float *dpar, *f, *bd_ax, *bd_bx, *bd_ay, *bd_by, *bd_az, *bd_bz;
  float  q;
  DFTI_DESCRIPTOR_HANDLE xhandle = 0;
  DFTI_DESCRIPTOR_HANDLE yhandle = 0;
  char *BCtype;

  
  dpar=(float*)malloc((5*(nx+ny)/2+9)*sizeof(float));
  f=(float*)malloc((nx+1)*(ny+1)*(nz+1)*sizeof(float));
  bd_ax=(float*)malloc((ny+1)*(nz+1)*sizeof(float));
  bd_bx=(float*)malloc((ny+1)*(nz+1)*sizeof(float));
  bd_ay=(float*)malloc((nx+1)*(nz+1)*sizeof(float));
  bd_by=(float*)malloc((nx+1)*(nz+1)*sizeof(float));
  bd_az=(float*)malloc((nx+1)*(ny+1)*sizeof(float));
  bd_bz=(float*)malloc((nx+1)*(ny+1)*sizeof(float));
  
  /* Defining the parallelepiped domain 0<x<1, 0<y<1, 0<z<1 for 3D Poisson Solver */
  ax=0.0E0;
  bx=1.0E0;
  ay=0.0E0;
  by=1.0E0;
  az=0.0E0;
  bz=1.0E0;

  /*******************************************************************************
  Setting the coefficient q to 0.
  Note that this is the way to use Helmholtz Solver to solve Poisson problem!
  *******************************************************************************/
  q=0.0E0;
  
  /* Computing the mesh size hx in x-direction */
  lx=bx-ax;
  hx=lx/nx;
  /* Computing the mesh size hy in y-direction */
  ly=by-ay;
  hy=ly/ny;
  /* Computing the mesh size hx in z-direction */
  lz=bz-az;
  hz=lz/nz;
  
  /*
  Here we are using the mesh sizes hx, hy, and hz computed before to compute
  the coordinates (xi,yi,zi) of the mesh points 
  */
  for(iz=0;iz<=nz;iz++)
  {
    for(iy=0;iy<=ny;iy++)
    {
      for(ix=0;ix<=nx;ix++)
      {
        xi=hx*ix/lx;
        yi=hy*iy/ly;
        zi=hz*iz/lz;
  
        f[ix+iy*(nx+1)+iz*(nx+1)*(ny+1)]=0;
      }
    }
  }
  BCtype = "DDDDDD";

  /* Setting the values of the boundary function G(x,y,z) that is equal to the TRUE solution
  in the mesh points laying on Dirichlet boundaries */
  for(iy=0;iy<=ny;iy++)
  {
    for(iz=0;iz<=nz;iz++)
    {
      bd_ax[iy+iz*(ny+1)]=1;
      bd_bx[iy+iz*(ny+1)]=2;
    }
  }
  /* Setting the values of the boundary function g(x,y,z) that is equal to the normal derivative
  of the TRUE solution in the mesh points laying on Neumann boundaries */
  for(ix=0;ix<=nx;ix++)
  {
    for(iz=0;iz<=nz;iz++)
    {
      bd_ay[ix+iz*(nx+1)]=3;
      bd_by[ix+iz*(nx+1)]=4;
    }
  }
  for(ix=0;ix<=nx;ix++)
  {
    for(iy=0;iy<=ny;iy++)
    {
      bd_az[ix+iy*(nx+1)]=5;
      bd_bz[ix+iy*(nx+1)]=6;
    }
  }
  
  /* Initializing ipar array to make it free from garbage */
  for(i=0;i<128;i++)
  {
    ipar[i]=0;
  }
  
  /* Initializing simple data structures of Poisson Library for 3D Poisson Solver */
  s_init_Helmholtz_3D(&ax, &bx, &ay, &by, &az, &bz, &nx, &ny, &nz, 
        BCtype, &q, ipar, dpar, &stat); CHKERRQ(stat);

  /* Initializing complex data structures of Poisson Library for 3D Poisson Solver
  NOTE: Right-hand side f may be altered after the Commit step. If you want to keep it,
  you should save it in another memory location! */
  s_commit_Helmholtz_3D(f, bd_ax, bd_bx, bd_ay, bd_by, bd_az, bd_bz, 
        &xhandle, &yhandle, ipar, dpar, &stat); CHKERRQ(stat);

  /* Computing the approximate solution of 3D Poisson problem
  NOTE: Boundary data stored in the arrays bd_ax, bd_bx, bd_ay, bd_by, bd_az, bd_bz should not be changed
          between the Commit step and the subsequent call to the Solver routine/*
          Otherwise the results may be wrong. */
  PetscBarrier(PETSC_NULL);
  PetscGetTime(&p1);
  s_Helmholtz_3D(f, bd_ax, bd_bx, bd_ay, bd_by, bd_az, bd_bz, 
        &xhandle, &yhandle, ipar, dpar, &stat); CHKERRQ(stat);
  PetscGetTime(&p2);
  /* Cleaning the memory used by xhandle and yhandle */
  free_Helmholtz_3D(&xhandle, &yhandle, ipar, &stat);
  /* Now we can use xhandle and yhandle to solve another 3D Poisson problem */
  
//  WriteVectorArray("res.Real64",(nx+1)*(ny+1)*(nz+1), f );
  
//  printf("\n\nTIME: %f\n\n", (p2-p1));
  *time = (p2 - p1);
  
  free(f);
  
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "SpectralMethod"
PetscInt EVENT_SpectralMethod;
PetscErrorCode SpectralMethod( int nn, PetscLogDouble *time )
{
  PetscErrorCode ierr;
  PetscLogDouble p1, p2;
  
  MKL_INT nx=nn, ny=nn, nz=nn;
  MKL_INT ix, iy, iz, i, stat;
  MKL_INT ipar[128];
  double ax, bx, ay, by, az, bz, lx, ly, lz, hx, hy, hz, xi, yi, zi, cx, cy, cz;
  double *dpar, *f, *bd_ax, *bd_bx, *bd_ay, *bd_by, *bd_az, *bd_bz;
  double q;
  DFTI_DESCRIPTOR_HANDLE xhandle = 0;
  DFTI_DESCRIPTOR_HANDLE yhandle = 0;
  char *BCtype;

  PetscFunctionBegin;
  PetscLogEventBegin(EVENT_SpectralMethod,0,0,0,0);

  dpar=(double*)malloc((5*(nx+ny)/2+9)*sizeof(double));
  f=(double*)malloc((nx+1)*(ny+1)*(nz+1)*sizeof(double));
  bd_ax=(double*)malloc((ny+1)*(nz+1)*sizeof(double));
  bd_bx=(double*)malloc((ny+1)*(nz+1)*sizeof(double));
  bd_ay=(double*)malloc((nx+1)*(nz+1)*sizeof(double));
  bd_by=(double*)malloc((nx+1)*(nz+1)*sizeof(double));
  bd_az=(double*)malloc((nx+1)*(ny+1)*sizeof(double));
  bd_bz=(double*)malloc((nx+1)*(ny+1)*sizeof(double));
  
  /* Defining the parallelepiped domain 0<x<1, 0<y<1, 0<z<1 for 3D Poisson Solver */
  ax=0.0E0;
  bx=1.0E0;
  ay=0.0E0;
  by=1.0E0;
  az=0.0E0;
  bz=1.0E0;

  /*******************************************************************************
  Setting the coefficient q to 0.
  Note that this is the way to use Helmholtz Solver to solve Poisson problem!
  *******************************************************************************/
  q=0.0E0;
  
  /* Computing the mesh size hx in x-direction */
  lx=bx-ax;
  hx=lx/nx;
  /* Computing the mesh size hy in y-direction */
  ly=by-ay;
  hy=ly/ny;
  /* Computing the mesh size hx in z-direction */
  lz=bz-az;
  hz=lz/nz;
  
  /*
  Here we are using the mesh sizes hx, hy, and hz computed before to compute
  the coordinates (xi,yi,zi) of the mesh points 
  */
  for(iz=0;iz<=nz;iz++)
  {
    for(iy=0;iy<=ny;iy++)
    {
      for(ix=0;ix<=nx;ix++)
      {
        xi=hx*ix/lx;
        yi=hy*iy/ly;
        zi=hz*iz/lz;
  
        f[ix+iy*(nx+1)+iz*(nx+1)*(ny+1)]=0;
      }
    }
  }
  BCtype = "DDDDDD";

  /* Setting the values of the boundary function G(x,y,z) that is equal to the TRUE solution
  in the mesh points laying on Dirichlet boundaries */
  for(iy=0;iy<=ny;iy++)
  {
    for(iz=0;iz<=nz;iz++)
    {
      bd_ax[iy+iz*(ny+1)]=1;
      bd_bx[iy+iz*(ny+1)]=2;
    }
  }
  /* Setting the values of the boundary function g(x,y,z) that is equal to the normal derivative
  of the TRUE solution in the mesh points laying on Neumann boundaries */
  for(ix=0;ix<=nx;ix++)
  {
    for(iz=0;iz<=nz;iz++)
    {
      bd_ay[ix+iz*(nx+1)]=3;
      bd_by[ix+iz*(nx+1)]=4;
    }
  }
  for(ix=0;ix<=nx;ix++)
  {
    for(iy=0;iy<=ny;iy++)
    {
      bd_az[ix+iy*(nx+1)]=5;
      bd_bz[ix+iy*(nx+1)]=6;
    }
  }
  
  /* Initializing ipar array to make it free from garbage */
  for(i=0;i<128;i++)
  {
    ipar[i]=0;
  }
  
  /* Initializing simple data structures of Poisson Library for 3D Poisson Solver */
  d_init_Helmholtz_3D(&ax, &bx, &ay, &by, &az, &bz, &nx, &ny, &nz, 
        BCtype, &q, ipar, dpar, &stat); CHKERRQ(stat);
ipar[22]=1;
  /* Initializing complex data structures of Poisson Library for 3D Poisson Solver
  NOTE: Right-hand side f may be altered after the Commit step. If you want to keep it,
  you should save it in another memory location! */
  d_commit_Helmholtz_3D(f, bd_ax, bd_bx, bd_ay, bd_by, bd_az, bd_bz, 
        &xhandle, &yhandle, ipar, dpar, &stat); CHKERRQ(stat);

  /* Computing the approximate solution of 3D Poisson problem
  NOTE: Boundary data stored in the arrays bd_ax, bd_bx, bd_ay, bd_by, bd_az, bd_bz should not be changed
          between the Commit step and the subsequent call to the Solver routine/*
          Otherwise the results may be wrong. */
  PetscBarrier(PETSC_NULL);
  PetscGetTime(&p1);
  d_Helmholtz_3D(f, bd_ax, bd_bx, bd_ay, bd_by, bd_az, bd_bz, 
        &xhandle, &yhandle, ipar, dpar, &stat); CHKERRQ(stat);
  PetscGetTime(&p2);
  /* Cleaning the memory used by xhandle and yhandle */
  free_Helmholtz_3D(&xhandle, &yhandle, ipar, &stat);
  /* Now we can use xhandle and yhandle to solve another 3D Poisson problem */
  
//  WriteVectorArray("res.Real64",(nx+1)*(ny+1)*(nz+1), f );
  
//  printf("\n\nTIME: %f\n\n", (p2-p1));
  *time = (p2 - p1);
  
  free(f);
  
  PetscLogEventEnd(EVENT_SpectralMethod,0,0,0,0);
  PetscFunctionReturn(0);
}




#undef __FUNCT__
#define __FUNCT__ "SpectralMethod2DTimer"
PetscInt EVENT_SpectralMethod2DTimer;
PetscErrorCode SpectralMethod2DTimer()
{
  printf("\n\nSpectralMethod2DTimer***************\n\n");
  PetscErrorCode ierr;
  int nn=14;
  MKL_INT nx=nn, ny=nn;           /* Note that the size of the transform nx must be even !!! */
  MKL_INT ix, iy, i, stat;
  MKL_INT ipar[128];
  double ax, bx, ay, by, lx, ly, hx, hy, xi, yi;
  double *dpar, *f, *bd_ax, *bd_bx, *bd_ay, *bd_by;
  double q;
  PetscLogDouble p1, p2;
  DFTI_DESCRIPTOR_HANDLE xhandle = 0;
  char *BCtype;
    
  PetscFunctionBegin;
  PetscLogEventBegin(EVENT_SpectralMethod2DTimer,0,0,0,0);

  dpar=(double*)malloc((13*nx/2+7)*sizeof(double));
  f=(double*)malloc((nx+1)*(ny+1)*sizeof(double));
  bd_ax=(double*)malloc((ny+1)*sizeof(double));
  bd_bx=(double*)malloc((ny+1)*sizeof(double));
  bd_ay=(double*)malloc((nx+1)*sizeof(double));
  bd_by=(double*)malloc((nx+1)*sizeof(double));

  /* Defining the rectangular domain 0<x<1, 0<y<1 for 2D Laplace Solver */
  ax=0.0E0;  bx=1.0E0;
  ay=0.0E0;  by=1.0E0;

  q=0.0E0;

  /* Computing the mesh size hx in x-direction */
  lx=bx-ax;
  hx=lx/nx;
  /* Computing the mesh size hy in y-direction */
  ly=by-ay;
  hy=ly/ny;

  for(iy=0;iy<=ny;iy++)
  {
    for(ix=0;ix<=nx;ix++)
    {
      xi=hx*ix/lx;
      yi=hy*iy/ly;
      f[ix+iy*(nx+1)]= -444.44444444444446;
    }
  }
  
  f[1] = 99555.55555555556;

  BCtype = "NNNN";
  
  for(iy=0;iy<=ny;iy++)
  {
    yi=hy*iy/ly;
    bd_ax[iy]=0;
    bd_bx[iy]=0;
  }
  
  for(ix=0;ix<=nx;ix++)
  {
    xi=hx*ix/lx;
    bd_ay[ix]=0;
    bd_by[ix]=0;
  }

  /* Initializing ipar array to make it free from garbage */
  for(i=0;i<128;i++)
  {
    ipar[i]=0;
  }

  d_init_Helmholtz_2D(&ax, &bx, &ay, &by, &nx, &ny, 
        BCtype, &q, ipar, dpar, &stat);CHKERRQ(stat);
  
  //Commit will alter RHS "f" and boundary values    
  d_commit_Helmholtz_2D(f, bd_ax, bd_bx, bd_ay, bd_by, 
        &xhandle, ipar, dpar, &stat); CHKERRQ(stat);
  

  PetscGetTime(&p1);
  d_Helmholtz_2D(f, bd_ax, bd_bx, bd_ay, bd_by, 
        &xhandle, ipar, dpar, &stat); //CHKERRQ(stat);
  PetscGetTime(&p2);
  
  free_Helmholtz_2D(&xhandle, ipar, &stat); CHKERRQ(stat);
  
  printf("TIME 2D: %f\n", (p2-p1) );

  WriteVectorArray("res.Real64",(nx+1)*(ny+1), f );
  
  PetscLogFlops(20*nn*nn*log(nn));
  PetscLogEventEnd(EVENT_SpectralMethod2DTimer,0,0,0,0);
  PetscFunctionReturn(0);
}

