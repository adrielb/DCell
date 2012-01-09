#include "ImmersedInterfaceMethod.h"
#include "IIM_private.h"

/* Surface quantities:
 *  -Updates normal, tangent, curvature.
 *  -Evaluates interfacial force function F
 */

#undef __FUNCT__
#define __FUNCT__ "IIMUpdateSurfaceQuantities"
PetscErrorCode IIMUpdateSurfaceQuantities( IIM iim, LevelSet ls )
{
  PetscErrorCode ierr;

  PetscFunctionBegin;

  if( ls->phi->is2D )
  {
    ierr = IIMUpdateIrregularNodeGrid_2D ( iim, ls ); CHKERRQ(ierr);
    ierr = IIMUpdateSurfaceQuantities_2D ( iim, ls ); CHKERRQ(ierr);
    ierr = IIMUpdateSurfaceDerivatives_2D( iim, ls ); CHKERRQ(ierr);
  } else {
    SETERRQ(PETSC_COMM_SELF,0,"NOT IMPLEMENTED\nIIMUpdateSurfaceQuantities_3D\nIIMUpdateSurfaceDerivatives_3D");
    /*
    ierr = IIMUpdateSurfaceQuantities_3D( iim, ls ); CHKERRQ(ierr);
    ierr = IIMUpdateSurfaceDerivatives_3D( iim, ls ); CHKERRQ(ierr);
    */
  }
  
  PetscFunctionReturn(0);
}

double LevelSetCurvature( Grid phi, Coor X, PetscReal dx )
{
  double px, py, px2, py2, pxx, pyy, pxy, k;
  double p[3][3];
  Coor x;
  Coor d = {dx,dx,0};
  int j,i;
  for (j = 0; j < 3; ++j) {
    for (i = 0; i < 3; ++i) {
      x.x = X.x - (i-1) * d.x;
      x.y = X.y - (j-1) * d.y;
      GridInterpolate( phi, x, &p[j][i] );
    }
  }
  i = j = 1;
  px = (p[j][i+1]-p[j][i-1]) / (2. * d.x);
  py = (p[j+1][i]-p[j-1][i]) / (2. * d.y);
  px2 = px * px;
  py2 = py * py;
  pxx = (p[j][i-1] - 2.*p[j][i] + p[j][i+1]) / (d.x*d.x);
  pyy = (p[j-1][i] - 2.*p[j][i] + p[j+1][i]) / (d.y*d.y);
  pxy = (p[j-1][i-1] + p[j+1][i+1] - p[j+1][i-1] - p[j-1][i+1]) / (4.*d.x*d.y);
  k   = (pxx*py2 - 2.*px*py*pxy + pyy*px2) / sqrt((px2+py2)*(px2+py2)*(px2+py2));
  return k;
}

#undef __FUNCT__
#define __FUNCT__ "IIMUpdateSurfaceQuantities_2D"
PetscErrorCode IIMUpdateSurfaceQuantities_2D( IIM iim, LevelSet ls )
{
  int i;
  IrregularNode *n;
  PetscReal h;
  PetscReal **phi;
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  ierr = PetscLogEventBegin(EVENT_IIMUpdateSurfaceQuantities,0,0,0,0); CHKERRQ(ierr);
  ierr = GridGet(ls->phi,&phi); CHKERRQ(ierr);
  iim->dh = ls->phi->d;
  for( i = 0; i < ArrayLength(ls->irregularNodes); i++ )
  {
    ierr = ArrayGet(ls->irregularNodes,i,&n); CHKERRQ(ierr);
    n->k = Bilinear2D( GridFunction2D_Curv, phi, iim->dh, n->X.x, n->X.y );
    n->k_nn = n->k;
    n->k_tt = 0;
    n->k_nt = 0;
    n->nx = Bilinear2D(GridFunction2D_DerivX, phi, iim->dh, n->X.x, n->X.y);
    n->ny = Bilinear2D(GridFunction2D_DerivY, phi, iim->dh, n->X.x, n->X.y);
    h = sqrt( n->nx*n->nx + n->ny*n->ny );
    n->nx = n->nx / h;
    n->ny = n->ny / h;
    n->sx = n->ny;
    n->sy =-n->nx;
    iim->F( n, iim->context );
  }
  ierr = PetscLogEventEnd(EVENT_IIMUpdateSurfaceQuantities,0,0,0,0); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "IIMUpdateSurfaceQuantities_3D"
PetscErrorCode IIMUpdateSurfaceQuantities_3D( IIM iim, LevelSet ls )
{
  PetscReal h;
  IrregularNode *n;
  Jump j; // Mangled, not a jump, but used for rotation functions
  PetscReal ***phi;
  Coor dh = ls->phi->d;
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  ierr = GridGet(ls->phi,&phi); CHKERRQ(ierr);
  int i;
  for( i = 0; i < ArrayLength(ls->irregularNodes); i++ )
  {
    ierr = ArrayGet(ls->irregularNodes,i,&n); CHKERRQ(ierr);
    
    n->nx = Bilinear3D(GridFunction3D_DerivX, phi, dh, n->X );
    n->ny = Bilinear3D(GridFunction3D_DerivY, phi, dh, n->X );
    n->nz = Bilinear3D(GridFunction3D_DerivZ, phi, dh, n->X );
    h = sqrt( n->nx*n->nx + n->ny*n->ny + n->nz*n->nz );
    n->nx /= h;
    n->ny /= h;
    n->nz /= h;
    LocalCoor3DTangential( n );
    
    j.x = n->nx * h;
    j.y = n->ny * h;
    j.z = n->nz * h;
    IIMGlobalToLocal_1st( n, &j); // inverted rotation! global (x,y,z) to local (e,n,t) coor sys! 
    j.xx = Bilinear3D(GridFunction3D_DerivXX, phi, dh, n->X );
    j.xy = Bilinear3D(GridFunction3D_DerivXY, phi, dh, n->X );
    j.xz = Bilinear3D(GridFunction3D_DerivXZ, phi, dh, n->X );
    j.yy = Bilinear3D(GridFunction3D_DerivYY, phi, dh, n->X );
    j.yz = Bilinear3D(GridFunction3D_DerivYZ, phi, dh, n->X );
    j.zz = Bilinear3D(GridFunction3D_DerivZZ, phi, dh, n->X );
    IIMGlobalToLocal_2nd( n, &j); 
    n->k_nn = -j.nn / j.e;
    n->k_tt = -j.tt / j.e;
    n->k_nt = -j.nt / j.e;
    n->k  = Bilinear3D(GridFunction3D_Curv,   phi, dh, n->X );
    //TODO: test curvatures values
    iim->F( n, iim->context );
  }
  //TODO: try to calculate # flops here
  
  PetscFunctionReturn(0);
}
