#include "mkl_poisson.h"
#include "IntelFFT.h"

#define CHKSTAT(stat) if( stat ) SETERRQ1(PETSC_ERR_LIB, "IntelFFT ERROR: %d", stat );
#define IPAR_SIZE 128

PetscErrorCode PCIntelFFT_Setup( void *ctx );
PetscErrorCode PCIntelFFT_Destroy( void *ctx );
PetscErrorCode PCIntelFFT_Apply(void *ctx, Vec xIn, Vec xOut);

struct _IntelFFT {
  MKL_INT *ipar;
  double ax, bx,
         ay, by,
         az, bz;
  double *bd_ax, *bd_bx, 
         *bd_ay, *bd_by, 
         *bd_az, *bd_bz;
  double *dpar, *f;
  double q;
  DFTI_DESCRIPTOR_HANDLE xhandle;
  DFTI_DESCRIPTOR_HANDLE yhandle;
  char *BCtype;
  
  int ov; // Overlap from PCASM
  IS is; // 1D global indexes from PCASM subdomain
  iCoor m; // global size of domain, number of grid points
  Coor d;  // Length btw grid points (in actual units)
  iCoor n; // local number of grid points in this fft
  iCoor s,e; // global coor of start and end index for fft 
  DALocalInfo g; // DA which contains grid block info  
};

#undef __FUNCT__
#define __FUNCT__ "PCIntelFFT_Setup"
PetscErrorCode PCIntelFFT_Setup( void *ctx )
{
  IntelFFT i = (IntelFFT) ctx;
  int *gxs = &i->g.xs, *gxm = &i->g.xm, *gmx = &i->g.mx, xe;
  int *s = &i->s.x, *e = &i->e.x;
  Coor d = i->d;
  iCoor n = i->n;
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  PetscInt ov = i->ov;
  int k;
  for( k = 0; k < 3; ++k)
  {
    int *n = &i->n.x;
    s[k] = gxs[k] == 0 ? s[k] : gxs[k] - i->ov;
    xe = gxs[k] + gxm[k];
    e[k] = xe == gmx[k] ? e[k] : xe + i->ov;
    n[k] = e[k] - s[k];
  }
  
  int db = sizeof(double),
      sx = db*(n.y+1)*(n.z+1),
      sy = db*(n.x+1)*(n.z+1),
      sz = db*(n.x+1)*(n.y+1),
      ss = (n.x+1)*(n.y+1)*(n.z+1);
  
  ierr = PetscNew(struct _IntelFFT,&i); CHKERRQ(ierr);
  ierr = PetscMalloc(128*sizeof(MKL_INT), &i->ipar); CHKERRQ(ierr);
  ierr = PetscMalloc(db*(5*(n.x+n.y)/2+9), &i->dpar); CHKERRQ(ierr);
  ierr = PetscMalloc(db*ss, &i->f); CHKERRQ(ierr);
  ierr = PetscMemzero(i->f,db*ss); CHKERRQ(ierr);
  
  ierr = PetscMalloc(sx,&i->bd_ax); CHKERRQ(ierr);
  ierr = PetscMalloc(sx,&i->bd_bx); CHKERRQ(ierr);
  ierr = PetscMemzero(i->bd_ax,sx); CHKERRQ(ierr);
  ierr = PetscMemzero(i->bd_bx,sx); CHKERRQ(ierr);
  
  ierr = PetscMalloc(sy,&i->bd_ay); CHKERRQ(ierr);
  ierr = PetscMalloc(sy,&i->bd_by); CHKERRQ(ierr);
  ierr = PetscMemzero(i->bd_ay,sy); CHKERRQ(ierr);
  ierr = PetscMemzero(i->bd_by,sy); CHKERRQ(ierr);
  
  ierr = PetscMalloc(sz,&i->bd_az); CHKERRQ(ierr);
  ierr = PetscMalloc(sz,&i->bd_bz); CHKERRQ(ierr);
  ierr = PetscMemzero(i->bd_az,sz); CHKERRQ(ierr);
  ierr = PetscMemzero(i->bd_bz,sz); CHKERRQ(ierr);
  
  i->ax = 0.0E0;   i->bx = n.x / d.x;
  i->ay = 0.0E0;   i->by = n.y / d.y;
  i->az = 0.0E0;   i->bz = n.z / d.z;
  i->xhandle = 0;  i->yhandle = 0;
  i->BCtype = "DDDDDD";
  i->q = 0.0E0;
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "PCIntelFFT_Destroy"
PetscErrorCode PCIntelFFT_Destroy( void *ctx )
{
  IntelFFT i = (IntelFFT) ctx;
  int stat;
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  
  if( i->xhandle && i->yhandle ) // Only if x/yhandle allocated
    free_Helmholtz_3D(&i->xhandle, &i->yhandle, i->ipar, &stat); CHKSTAT(stat);
  ierr = PetscFree(i->bd_ax); CHKERRQ(ierr);
  ierr = PetscFree(i->bd_bx); CHKERRQ(ierr);
  ierr = PetscFree(i->bd_ay); CHKERRQ(ierr);
  ierr = PetscFree(i->bd_by); CHKERRQ(ierr);
  ierr = PetscFree(i->bd_az); CHKERRQ(ierr);
  ierr = PetscFree(i->bd_bz); CHKERRQ(ierr);
  ierr = PetscFree(i->dpar); CHKERRQ(ierr);
  ierr = PetscFree(i->ipar); CHKERRQ(ierr);
  ierr = PetscFree(i->f); CHKERRQ(ierr);
  ierr = PetscFree(i); CHKERRQ(ierr);
  
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "IntelFFTSolve"
PetscInt EVENT_IntelFFTSolve;
PetscErrorCode IntelFFTSolve( IntelFFT i )
{
  int stat;
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  PetscLogEventBegin(EVENT_IntelFFTSolve,0,0,0,0);
//  PetscLogEventRegister(&EVENT_IntelFFTSolve,"IntelFFTSolve", 0);
  
  d_init_Helmholtz_3D(
      &i->ax, &i->bx, 
      &i->ay, &i->by, 
      &i->az, &i->bz, 
      &i->n.x, &i->n.y, &i->n.z, 
      i->BCtype, &i->q, i->ipar, i->dpar, &stat); CHKSTAT(stat);
  
  d_commit_Helmholtz_3D( i->f, 
      i->bd_ax, i->bd_bx, 
      i->bd_ay, i->bd_by, 
      i->bd_az, i->bd_bz, 
      &i->xhandle, &i->yhandle, 
      i->ipar, i->dpar, &stat); CHKSTAT(stat);
  
  d_Helmholtz_3D( i->f, 
      i->bd_ax, i->bd_bx, 
      i->bd_ay, i->bd_by, 
      i->bd_az, i->bd_bz, 
      &i->xhandle, &i->yhandle, 
      i->ipar, i->dpar, &stat); CHKSTAT(stat);

  PetscLogEventEnd(EVENT_IntelFFTSolve,0,0,0,0);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PCIntelFFTSetPC"
PetscErrorCode PCIntelFFTSetPC( PC pc, DALocalInfo g, IS is, iCoor s, iCoor e )
{
  IntelFFT fft;
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  ierr = PetscMalloc(sizeof(struct _IntelFFT),&fft); CHKERRQ(ierr);
  fft->g = g;
  fft->is = is;
  fft->s = s;
  fft->e = e;
  
  ierr = PCSetType(pc,PCSHELL); CHKERRQ(ierr);
  ierr = PCShellSetContext(pc,fft); CHKERRQ(ierr);
  ierr = PCShellSetApply(pc,PCIntelFFT_Apply); CHKERRQ(ierr);
  ierr = PCShellSetSetUp(pc,PCIntelFFT_Setup); CHKERRQ(ierr);
  ierr = PCShellSetDestroy(pc,PCIntelFFT_Destroy); CHKERRQ(ierr);
  
  PetscFunctionReturn(0);
}



#undef __FUNCT__
#define __FUNCT__ "PCIntelFFT_Apply"
PetscErrorCode PCIntelFFT_Apply(void *ctx, Vec xIn, Vec xOut)
{
  PetscInt len, g, ix, iy, iz;
  PetscReal *xin, *xout;
  const PetscInt *is;
  IntelFFT fft = (IntelFFT) ctx;
  iCoor s = fft->s, e = fft->e, m = fft->m, n = fft->n;
  PetscInt gx, gy, gz;
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
      
  //Scatter values from xin to F
  ierr = ISGetSize(fft->is,&len); CHKERRQ(ierr);
  ierr = ISGetIndices(fft->is,&is); CHKERRQ(ierr);
  ierr = VecGetArray(xIn, &xin); CHKERRQ(ierr);
  int q;
  for( q = 0; q < len; ++q) // q: local 1D index 
  {
    g = is[q]; // g: global 1D index
    // convert 'g' to global coor (gx, gy, gz)
    gz =   g              / ( m.x*m.y );
    gy = ( g - gz*m.x*m.y ) /   m.x;
    gx =   g - gz*m.x*m.y - gy*m.x;
    
    //filter coor outside fft box
    if( s.x < gx && gx < e.x && 
        s.y < gy && gy < e.y &&
        s.z < gz && gz < e.z )
    {
      // (xsf, ysf, zsf) global coor of fft origin
      // (ix, iy, iz) local fft coor
      ix = gx - s.x;
      iy = gy - s.y;
      iz = gz - s.z;
      
      fft->f[ix+iy*(n.x+1)+iz*(n.x+1)*(n.y+1)] = xin[q];
    }
  }
  ierr = VecRestoreArray(xIn, &xin); CHKERRQ(ierr);
  
  //Solve FFT(F)
  ierr = IntelFFTSolve( fft ); CHKERRQ(ierr);
  
  //Scatter values from F to xout
  ierr = VecGetArray(xOut, &xout); CHKERRQ(ierr);
  for( q = 0; q < len; ++q) // q: local 1D index 
  {
    g = is[q]; // g: global 1D index
    // convert 'g' to global coor (gx, gy, gz)
    gz =   g              / ( m.x*m.y );
    gy = ( g - gz*m.x*m.y ) /   m.x;
    gx =   g - gz*m.x*m.y - gy*m.x;
    
    //filter coor outside fft box
    if( s.x < gx && gx < e.x && 
        s.y < gy && gy < e.y &&
        s.z < gz && gz < e.z )
    {
      xout[q] =fft->f[ix+iy*(n.x+1)+iz*(n.x+1)*(n.y+1)];
    }
  }
  ierr = VecRestoreArray(xOut, &xout); CHKERRQ(ierr);
  ierr = ISRestoreIndices(fft->is,&is); CHKERRQ(ierr);
  
  PetscFunctionReturn(0);
}