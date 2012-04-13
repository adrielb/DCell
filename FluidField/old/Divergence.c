#include "FluidField.h"

#undef __FUNCT__
#define __FUNCT__ "FluidFieldDivergence"
PetscLogEvent EVENT_FluidFieldDivergence;
PetscErrorCode FluidFieldDivergence(FluidField f)
{
  DALocalInfo info;
  Vec u,v,w; // Local, ghosted vectors
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  PetscLogEventBegin(EVENT_FluidFieldDivergence,0,0,0,0);
//  PetscLogEventRegister(&EVENT_FluidFieldDivergence,"FluidFieldDivergence", 0);
  ierr = DAGetLocalInfo(f->da,&info); CHKERRQ(ierr);
  ierr = DAGetLocalVector(f->da,&u); CHKERRQ(ierr);
  ierr = DAGetLocalVector(f->da,&v); CHKERRQ(ierr);
  ierr = DAGlobalToLocalBegin(f->da,f->u,INSERT_VALUES,u); CHKERRQ(ierr);
  ierr = DAGlobalToLocalEnd(  f->da,f->u,INSERT_VALUES,u); CHKERRQ(ierr);
  ierr = DAGlobalToLocalBegin(f->da,f->v,INSERT_VALUES,v); CHKERRQ(ierr);
  //TODO: interleaving work with communication here a possible source of optimization?
  ierr = DAGlobalToLocalEnd(  f->da,f->v,INSERT_VALUES,v); CHKERRQ(ierr);
  
  if( f->is3D )
  {
    ierr = DAGetLocalVector(f->da,&w); CHKERRQ(ierr);
    ierr = DAGlobalToLocalBegin(f->da,f->w,INSERT_VALUES,w); CHKERRQ(ierr);
    ierr = DAGlobalToLocalEnd(  f->da,f->w,INSERT_VALUES,w); CHKERRQ(ierr);
    FluidFieldDivergence_3D(info, f->d, u, v, w, f->div);
    ierr = DARestoreLocalVector(f->da,&w); CHKERRQ(ierr);
  } else {
    FluidFieldDivergence_2D(info, f->d, u, v, f->div);
  }
  ierr = DARestoreLocalVector(f->da,&u); CHKERRQ(ierr);
  ierr = DARestoreLocalVector(f->da,&v); CHKERRQ(ierr);
  
  PetscLogEventEnd(EVENT_FluidFieldDivergence,0,0,0,0);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "FluidFieldDivergence_2D"
PetscErrorCode FluidFieldDivergence_2D( DALocalInfo info, Coor d, 
    Vec US, Vec VS, Vec DIV )
{
  int i,j;
  PetscReal xe, xs, ye, ys;
  PetscReal **u, **v, **div;
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
  
  ierr = DAVecGetArray(info.da,US,&u); CHKERRQ(ierr);
  ierr = DAVecGetArray(info.da,VS,&v); CHKERRQ(ierr);
  ierr = DAVecGetArray(info.da,DIV,&div); CHKERRQ(ierr);
  for (j = ys; j < ye; ++j)
  {
    for (i = xs; i < xe; ++i)
    {
      div[j][i] = (u[j][i+1] - u[j][i]) / d.x +
                  (v[j+1][i] - v[j][i]) / d.y;
    }
  }
  ierr = DAVecRestoreArray(info.da,US,&u); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(info.da,VS,&v); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(info.da,DIV,&div); CHKERRQ(ierr);
  PetscLogFlops(5*(xe-xs)*(ye-ys));
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "FluidFieldDivergence_3D"
PetscErrorCode FluidFieldDivergence_3D( DALocalInfo info,
    Coor d, Vec US, Vec VS, Vec WS, Vec DIV )
{
  int i,j,k;
  PetscReal xe, xs, ye, ys, ze, zs;
  PetscReal ***u, ***v, ***w, ***div;
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
  ierr = DAVecGetArray(info.da,US,&u); CHKERRQ(ierr);
  ierr = DAVecGetArray(info.da,VS,&v); CHKERRQ(ierr);
  ierr = DAVecGetArray(info.da,WS,&w); CHKERRQ(ierr);
  ierr = DAVecGetArray(info.da,DIV,&div); CHKERRQ(ierr);
  for( k = zs; k < ze; ++k)
  {
    for (j = ys; j < ye; ++j)
    {
      for (i = xs; i < xe; ++i)
      {
        div[k][j][i] = (u[k][j][i+1] - u[k][j][i]) / d.x +
                       (v[k][j+1][i] - v[k][j][i]) / d.y +
                       (w[k+1][j][i] - w[k][j][i]) / d.z;
      }
    }
  }
  ierr = DAVecRestoreArray(info.da,US,&u); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(info.da,VS,&v); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(info.da,WS,&w); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(info.da,DIV,&div); CHKERRQ(ierr);
  PetscLogFlops(8*(xe-xs)*(ye-ys)*(ze-zs));
  PetscFunctionReturn(0);
}