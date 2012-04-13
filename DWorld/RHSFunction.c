#include "DWorld.h"

#undef __FUNCT__
#define __FUNCT__ "DWorldRHSFunction2D"
PetscInt EVENT_DWorldRHSFunction2D;
PetscErrorCode DWorldRHSFunction2D( TS ts,PetscReal t,Vec globalin,Vec globalout,void *ptr)
{
  DWorld ctx = (DWorld)ptr;
  DA da = ctx->da;
  Reaction rxn = ctx->rxnWorld;
  PetscInt xe,ye, xl, xr, yl, yr;
  PetscReal D;
  PetscReal ***c, ***chem;
  PetscReal dx2 = ctx->d.x * ctx->d.x, 
            dy2 = ctx->d.y * ctx->d.y; 
  int i, j, cc;
  Vec local;
  DALocalInfo l;
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  PetscLogEventBegin(EVENT_DWorldRHSFunction2D,0,0,0,0);
//  PetscLogEventRegister(&EVENT_DWorldRHSFunction2D,"DWorldRHSFunction2D", 0);
  ierr = DAGetLocalVector(da,&local); CHKERRQ(ierr);
  ierr = DAGlobalToLocalBegin(da,globalin,INSERT_VALUES,local);CHKERRQ(ierr);
  ierr = DAGlobalToLocalEnd  (da,globalin,INSERT_VALUES,local);CHKERRQ(ierr);
  ierr = DAGetLocalInfo(da,&l); CHKERRQ(ierr);
 
  xe = l.xs + l.xm;  // End index in local domain
  ye = l.ys + l.ym;

  ierr = DAVecGetArrayDOF(da,local,&c);CHKERRQ(ierr);
  ierr = DAVecGetArrayDOF(da,globalout,&chem);CHKERRQ(ierr);
  
  for(j=l.ys; j<ye; j++) 
  {
    for(i=l.xs; i<xe; i++)
    {
      // Enforce no-flux boundary condition
      xl = i == 0    ? i : i - 1;
      xr = i == l.mx-1 ? i : i + 1;
      yl = j == 0    ? j : j - 1;
      yr = j == l.my-1 ? j : j + 1; 
      
      ReactionUpdateFunction( rxn, c[j][i] );
      
      for (cc = 0; cc < l.dof; ++cc)
      {
        D = rxn->D[cc];
        
        chem[j][i][cc] = D / dy2 * (c[yl][i][cc] + c[yr][i][cc]) +
                         D / dx2 * (c[j][xl][cc] + c[j][xr][cc]) -
                         D * 2 * ( 1/dx2 + 1/dy2 ) * c[j][i][cc] +
                         rxn->F[cc];
      }
    }
  }

  ierr = DAVecRestoreArrayDOF(da,local,&c);CHKERRQ(ierr);
  ierr = DAVecRestoreArrayDOF(da,globalout,&chem);CHKERRQ(ierr);
  ierr = DARestoreLocalVector(da,&local); CHKERRQ(ierr);
  
  PetscLogEventEnd(EVENT_DWorldRHSFunction2D,0,0,0,0);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "RHSFunction"
PetscInt EVENT_RHSFunction;
PetscErrorCode RHSFunction(TS ts,PetscReal t,Vec vec_c,Vec vec_chem,void *ptr)
{
  DWorld ctx = (DWorld)ptr;
  DA da = ctx->da;
  Reaction rxn = ctx->rxnWorld;
  PetscInt xm, ym, zm, xs, ys, zs, mx, my, mz, dof, xe, ye, ze, xl, xr, yl, yr, zl, zr;
  PetscReal dx, dy, dz, dx2, dy2, dz2;
  PetscReal D;
  PetscReal ****c, ****chem;
  int i, j, k, cc;
  Vec local;
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  PetscLogEventBegin(EVENT_RHSFunction,0,0,0,0);
  ierr = DAGetLocalVector(da,&local); CHKERRQ(ierr);
  ierr = DAGlobalToLocalBegin(da,vec_c,INSERT_VALUES,local);CHKERRQ(ierr);
  ierr = DAGlobalToLocalEnd  (da,vec_c,INSERT_VALUES,local);CHKERRQ(ierr);
  ierr = DAGetInfo(da,0,&mx,&my,&mz,0,0,0,&dof,0,0,0);CHKERRQ(ierr);
  dx = ctx->len.x / (mx-1);
  dx2= dx*dx;
  dy = ctx->len.y / (my-1);
  dy2= dy*dy;
  dz = ctx->len.z / (mz-1);
  dz2= dz*dz;
  ierr = DAGetCorners(da,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
  xe = xs + xm;  // End index in local domain
  ye = ys + ym;
  ze = zs + zm;
  ierr = DAVecGetArrayDOF(da,local,&c);CHKERRQ(ierr);
  ierr = DAVecGetArrayDOF(da,vec_chem,&chem);CHKERRQ(ierr);
  
  for(k=zs; k<ze; k++)
  {
    for(j=ys; j<ye; j++) 
    {
      for(i=xs; i<xe; i++)
      {
        // Enforce no-flux boundary condition
        xl = i == 0    ? i : i - 1;
        xr = i == mx-1 ? i : i + 1;
        yl = j == 0    ? j : j - 1;
        yr = j == my-1 ? j : j + 1; 
        zl = k == 0    ? k : k - 1; 
        zr = k == mz-1 ? k : k + 1;
        
        /*
        for s in cells // each dcell
        {
          if( (i,j,k) in cell 's' )
            isExtracellular = false;
        }
        
        if( isExtracelluar )
        {*/
          ReactionUpdateFunction( rxn, c[k][j][i] );
          
          for (cc = 0; cc < dof; ++cc)
          {
            D = rxn->D[cc];
            
            chem[k][j][i][cc] = D / dz2 * (c[zl][j][i][cc] + c[zr][j][i][cc]) +
                                D / dy2 * (c[k][yl][i][cc] + c[k][yr][i][cc]) +
                                D / dx2 * (c[k][j][xl][cc] + c[k][j][xr][cc]) -
                                D * 2 * ( 1/dx2 + 1/dy2 + 1/dz2) * c[k][j][i][cc] +
                                rxn->F[cc];
          }
        
      }
    }
  }
  /*
  for s in cells // each dcell
  {
    DCellRHSFunction( cell[s], da, vec_in, vec_out );
    for( (i,j,k) inside cell[s] && (i,j,k) inside local domain )
    {
      update vec_out
    }
  }*/

  ierr = DAVecRestoreArrayDOF(da,local,&c);CHKERRQ(ierr);
  ierr = DAVecRestoreArrayDOF(da,vec_chem,&chem);CHKERRQ(ierr);
  ierr = DARestoreLocalVector(da,&local); CHKERRQ(ierr);
  PetscLogEventEnd(EVENT_RHSFunction,0,0,0,0);
  PetscFunctionReturn(0);
}