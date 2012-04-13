#include "FluidField.h"

#undef __FUNCT__
#define __FUNCT__ "FluidFieldGradient"
PetscInt EVENT_FluidFieldGradient;
PetscErrorCode FluidFieldGradient( FluidField f )
{
  Vec p_local, s_local;
  DALocalInfo info;
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  PetscLogEventBegin(EVENT_FluidFieldGradient,0,0,0,0);
//  PetscLogEventRegister(&EVENT_FluidFieldGradient,"FluidFieldGradient", 0);
  
  ierr = DAGetLocalVector(f->da,&p_local); CHKERRQ(ierr);
  ierr = DAGetLocalVector(f->da,&s_local); CHKERRQ(ierr);
  ierr = DAGlobalToLocalBegin(f->da,f->p,INSERT_VALUES,p_local); CHKERRQ(ierr);
  ierr = DAGlobalToLocalEnd(  f->da,f->p,INSERT_VALUES,p_local); CHKERRQ(ierr);
  ierr = DAGlobalToLocalBegin(f->da,f->source,INSERT_VALUES,s_local); CHKERRQ(ierr);
  ierr = DAGlobalToLocalEnd(  f->da,f->source,INSERT_VALUES,s_local); CHKERRQ(ierr);
  ierr = DAGetLocalInfo(f->da,&info); CHKERRQ(ierr);
  if( f->is3D )
  {
    FluidFieldGradient_3D( info, f->d, f->mu, f->rho, p_local, s_local, f->px, f->py, f->pz );
  } else {
    FluidFieldGradient_2D( info, f->d, f->mu, f->rho, p_local, s_local, f->px, f->py );
  }
  ierr = DARestoreLocalVector(f->da,&p_local); CHKERRQ(ierr);
  ierr = DARestoreLocalVector(f->da,&s_local); CHKERRQ(ierr);
  
  PetscLogEventEnd(EVENT_FluidFieldGradient,0,0,0,0);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "FluidFieldGradient_2D"
PetscErrorCode FluidFieldGradient_2D(DALocalInfo info, Coor d, PetscReal mu, PetscReal rho, 
    Vec P, Vec S, Vec PX, Vec PY )
{
  int i, j;
  PetscReal xe, xs, ye, ys;
  PetscReal **s, **p, **px, **py;
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  xs = info.xs;
  xs = xs == 0 ? 2 : xs;
  xe = info.xs+info.xm;
  xe = xe == info.mx ? info.mx-3 : xe;
  ys = info.ys;
  ys = ys == 0 ? 2 : ys;
  ye = info.ys+info.ym;
  ye = ye == info.my ? info.my-3 : ye;
  
  DAVecGetArray(info.da,P,&p);
  DAVecGetArray(info.da,S,&s);
  DAVecGetArray(info.da,PX,&px);
  DAVecGetArray(info.da,PY,&py);
  for( j = ys; j < ye; ++j)
  {
    for( i = xs; i < xe; ++i)
    {
      px[j][i] = ( p[j][i] - p[j][i-1] ) / (d.x * mu) - ( s[j][i] - s[j][i-1] ) / (d.x * 3 * rho);
      py[j][i] = ( p[j][i] - p[j-1][i] ) / (d.y * mu) - ( s[j][i] - s[j-1][i] ) / (d.y * 3 * rho);
      if( j==1 || j==info.my-2 ) py[j][i] = 0;
      if( i==1 || i==info.mx-2 ) px[j][i] = 0;
    }
  }
  DAVecRestoreArray(info.da,P,&p);
  DAVecRestoreArray(info.da,S,&s);
  DAVecRestoreArray(info.da,PX,&px);
  DAVecRestoreArray(info.da,PY,&py); 
  PetscLogFlops(5*2*(xe-xs)*(ye-ys)); // 5 flops per dim = 10
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "FluidFieldGradient_3D"
PetscErrorCode FluidFieldGradient_3D(DALocalInfo info, Coor d, PetscReal mu, PetscReal rho,
    Vec P, Vec S, Vec PX, Vec PY, Vec PZ )
{
  int i, j, k;
  PetscReal xe, xs, ye, ys, ze, zs;
  PetscReal ***s, ***p, ***px, ***py, ***pz;
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  
  xs = info.xs;
  xs = xs == 0 ? 2 : xs;
  xe = info.xs+info.xm;
  xe = xe == info.mx ? info.mx-3 : xe;
  ys = info.ys;
  ys = ys == 0 ? 2 : ys;
  ye = info.ys+info.ym;
  ye = ye == info.my ? info.my-3 : ye;
  zs = info.zs;
  zs = zs == 0 ? 2 : zs;
  ze = info.zs+info.zm;
  ze = ze == info.mz ? info.mz-3 : ze;
  
  DAVecGetArray(info.da,S,&s);
  DAVecGetArray(info.da,P,&p);
  DAVecGetArray(info.da,PX,&px);
  DAVecGetArray(info.da,PY,&py);
  DAVecGetArray(info.da,PZ,&pz);
  for( k = zs; k < ze; ++k)
  {
    for( j = ys; j < ye; ++j)
    {
      for( i = xs; i < xe; ++i)
      {
        px[k][j][i] = ( p[k][j][i] - p[k][j][i-1] ) / (d.x * mu) - ( s[k][j][i] - s[k][j][i-1] ) / (d.x * 3 * rho);
        py[k][j][i] = ( p[k][j][i] - p[k][j-1][i] ) / (d.y * mu) - ( s[k][j][i] - s[k][j-1][i] ) / (d.y * 3 * rho);
        pz[k][j][i] = ( p[k][j][i] - p[k-1][j][i] ) / (d.z * mu) - ( s[k][j][i] - s[k-1][j][i] ) / (d.z * 3 * rho);
        if( i==1 || i==info.mx-2 ) px[k][j][i] = 0;
        if( j==1 || j==info.my-2 ) py[k][j][i] = 0;
        if( k==1 || k==info.mz-2 ) pz[k][j][i] = 0;
      }
    }
  }
  DAVecRestoreArray(info.da,S,&s);
  DAVecRestoreArray(info.da,P,&p);
  DAVecRestoreArray(info.da,PX,&px);
  DAVecRestoreArray(info.da,PY,&py);
  DAVecRestoreArray(info.da,PZ,&pz);
  PetscLogFlops(5*3*(xe-xs)*(ye-ys)*(ze-zs)); // 5 flops per dim = 15
  PetscFunctionReturn(0);
}