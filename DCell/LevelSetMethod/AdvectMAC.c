#include "LevelSetMethod.h"

  /*           .-----vij+1---.
   *           |      |      |
   * pij-1----uij----pij----uij+1----pij+1
   *           |      |      |
   *           .-----vij-----.
   */
//TODO: advection of tube instead of the entire domain
//TODO: extension velocity: u'[t] + F | Del( u ) |
/* Conservative form
 * u'[t] + ( a[x] . u )_x == 0
 */

#undef __FUNCT__
#define __FUNCT__ "LevelSetAdvectMAC"
PetscInt EVENT_LevelSetAdvectMAC;
PetscErrorCode LevelSetAdvectMAC( LevelSet ls, Coor d, PetscReal dt, 
    PetscReal **u, PetscReal **v, LevelSet new )
{
  int i, j, xs, xe, ys, ye;
  PetscReal **phi = ls->g->v2, **fi = new->g->v2;
  PetscReal p, am, an, bm, bn, du, dv;
  iCoor pos = ls->g->p;
  iCoor num = ls->g->n;
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  ierr = PetscLogEventBegin(EVENT_LevelSetAdvectMAC,0,0,0,0); CHKERRQ(ierr);
//  PetscLogEventRegister(&EVENT_LevelSet2DAdvectMAC,"LevelSet2DAdvectMAC", 0);
  
  xs = pos.x < 0 ? 0 : pos.x;
  xe = pos.x + num.x - 1;
//  xe = xe > ? : xe; TODO: xs, xe in velocity field may not correspond to level set indexing
  
  ys = pos.y < 0 ? 0 : pos.y;
  ye = pos.y + num.y - 1;
  
  
  for (j = ys; j < ye; ++j)
  {
    for (i = xs; i < xe; ++i)
    {
      // [m----n]
      p = phi[j][i];
      am = i == 0    ? 2*p - phi[j][i+1] : phi[j][i-1];
      an = i == xe-1 ? 2*p - phi[j][i-1] : phi[j][i+1];
      bm = j == 0    ? 2*p - phi[j+1][i] : phi[j-1][i];
      bn = j == ye-1 ? 2*p - phi[j-1][i] : phi[j+1][i]; 

      du = u[j][i]   > 0. ? am * u[j][i]   :  p * u[j][i];
      dv = v[j][i]   > 0. ? bm * v[j][i]   :  p * v[j][i];
      du-= u[j][i+1] > 0. ?  p * u[j][i+1] : an * u[j][i+1];
      dv-= v[j+1][i] > 0. ?  p * v[j+1][i] : bn * v[j+1][i];

      fi[j][i] = phi[j][i] + dt * (du / d.x + dv / d.y);
    }
  }
  ierr = PetscLogEventEnd(EVENT_LevelSetAdvectMAC,0,0,0,0); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "LevelSetAdvectMAC_3D"
PetscErrorCode LevelSetAdvectMAC_3D(LevelSet ls, Coor d, PetscReal dt, 
    PetscReal ***u, PetscReal ***v, PetscReal ***w, LevelSet new )
{
  int i, j, k, xs, xe, ys, ye, zs, ze;
  PetscReal ***phi = ls->g->v3, ***fi = new->g->v3;
  PetscReal p, am, an, bm, bn, cm, cn, du, dv, dw;
  iCoor pos = ls->g->p;
  iCoor num = ls->g->n;
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  xs = pos.x < 0 ? 0 : pos.x;
  xe = pos.x + num.x - 1;
//  xe = xe > ? : xe; TODO: xs, xe in velocity field may not correspond to level set indexing
  ys = pos.y < 0 ? 0 : pos.y;
  ye = pos.y + num.y - 1;
  zs = pos.z < 0 ? 0 : pos.z;
  ze = pos.z + num.z - 1;
    
  for( k = zs; k < ze; ++k)
  {
    for (j = ys; j < ye; ++j)
    {
      for (i = xs; i < xe; ++i)
      {
        // [m----n]
        p = phi[k][j][i];
        am = i == 0    ? 2*p - phi[k][j][i+1] : phi[k][j][i-1];
        an = i == xe-1 ? 2*p - phi[k][j][i-1] : phi[k][j][i+1];
        bm = j == 0    ? 2*p - phi[k][j+1][i] : phi[k][j-1][i];
        bn = j == ye-1 ? 2*p - phi[k][j-1][i] : phi[k][j+1][i];
        cm = k == 0    ? 2*p - phi[k+1][j][i] : phi[k-1][j][i];
        cn = k == ze-1 ? 2*p - phi[k-1][j][i] : phi[k+1][j][i];

        du = u[k][j][i]   > 0. ? am * u[k][j][i]   :  p * u[k][j][i];
        dv = v[k][j][i]   > 0. ? bm * v[k][j][i]   :  p * v[k][j][i];
        dw = w[k][j][i]   > 0. ? cm * w[k][j][i]   :  p * w[k][j][i];
        du-= u[k][j][i+1] > 0. ?  p * u[k][j][i+1] : an * u[k][j][i+1];
        dv-= v[k][j+1][i] > 0. ?  p * v[k][j+1][i] : bn * v[k][j+1][i];
        dw-= w[k+1][j][i] > 0. ?  p * w[k+1][j][i] : cn * w[k+1][j][i];
        
        fi[k][j][i] = phi[k][j][i] + dt * (du / d.x + dv / d.y + dw / d.z);
      }
    }
  }
  PetscFunctionReturn(0);
}
  
  

/*
#undef __FUNCT__
#define __FUNCT__ "LevelSet2DAdvectMAC"
PetscInt EVENT_LevelSet2DAdvectMAC;
PetscErrorCode LevelSet2DAdvectMAC( LevelSet2D ls, Coor dh, PetscReal dt, 
    PetscReal **u, PetscReal **v, LevelSet2D new )
{
  int i, j;
  PetscReal **phi = ls->g2d->v2, **fi = new->g2d->v2;
  PetscReal dpx, dpy;
  PetscReal dx = 2*dh.x, dy = 2*dh.y;
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  PetscLogEventBegin(EVENT_LevelSet2DAdvectMAC,0,0,0,0);
//  PetscLogEventRegister(&EVENT_LevelSet2DAdvectMAC,"LevelSet2DAdvectMAC", 0);
  
  for (i = 1; i < ls->g2d->d1 - 1; ++i)
  {
    for (j = 1; j < ls->g2d->d2 - 1; ++j)
    {
      dpx = u[i][j] * ( phi[i-1][j]+phi[i][j] ) - u[i+1][j] * ( phi[i][j]+phi[i+1][j] );
      dpy = v[i][j] * ( phi[i][j-1]+phi[i][j] ) - v[i][j+1] * ( phi[i][j]+phi[i][j+1] ); 
      fi[i][j] = phi[i][j] + dt * (dpx / dx + dpy / dy);
    }
  }
  
  //TODO: figure out boundary conditions for 2nd MAC advection
  //TODO: try out higher order time integration beside explicit euler?
  //TODO: without reinitialization, 2nd excessive oscillations are present
  
  PetscLogEventEnd(EVENT_LevelSet2DAdvectMAC,0,0,0,0);
  PetscFunctionReturn(0);
}
*/
