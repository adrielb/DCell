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
  ierr = IIMUpdateIrregularNodes( iim, ls ); CHKERRQ(ierr);

  if( ls->phi->is2D )
  {
    ierr = IIMUpdateSurfaceQuantities_2D ( iim, ls ); CHKERRQ(ierr);
    ierr = IIMUpdateSurfaceDerivatives_2D( iim, ls ); CHKERRQ(ierr);
  } else {
    ierr = IIMUpdateSurfaceQuantities_3D ( iim, ls ); CHKERRQ(ierr);
    ierr = IIMUpdateSurfaceDerivatives_3D( iim, ls ); CHKERRQ(ierr);
  }
  
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "IIMUpdateSurfaceQuantities_2D"
PetscErrorCode IIMUpdateSurfaceQuantities_2D( IIM iim, LevelSet ls )
{
  int i;
  IIMIrregularNode *n;
  Coor nv;
  const int len = ArrayLength( iim->iimIrregularNodes );
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  ierr = PetscLogEventBegin(EVENT_IIMUpdateSurfaceQuantities,0,0,0,0); CHKERRQ(ierr);
  for( i = 0; i < len; i++ )
  {
    ierr = ArrayGet( iim->iimIrregularNodes ,i,&n); CHKERRQ(ierr);
    n->k = GridBilinear( ls->phi, GridFunction2D_Curv, n->X );
    n->k_nn = n->k;
    n->k_tt = 0;
    n->k_nt = 0;
    ierr = LevelSetNormalDirection( ls, n->X, &nv ); CHKERRQ(ierr);
    n->nx =  nv.x;
    n->ny =  nv.y;
    n->sx =  nv.y;
    n->sy = -nv.x;
    iim->F( n, iim->context );
  }
  ierr = PetscLogEventEnd(EVENT_IIMUpdateSurfaceQuantities,0,0,0,0); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "IIMUpdateSurfaceQuantities_3D"
PetscErrorCode IIMUpdateSurfaceQuantities_3D( IIM iim, LevelSet ls )
{
  int i;
  const int len = ArrayLength( iim->iimIrregularNodes );
  IIMIrregularNode *n;
  Jump j; // Mangled, not a jump, but used for rotation functions
  Coor nv;
  Grid phi = ls->phi;
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  ierr = PetscLogEventBegin(EVENT_IIMUpdateSurfaceQuantities,0,0,0,0); CHKERRQ(ierr);
  for( i = 0; i < len; i++ )
  {
    ierr = ArrayGet( iim->iimIrregularNodes, i,&n); CHKERRQ(ierr);
    
    // Surface Normal
    ierr = LevelSetNormalDirection( ls, n->X, &nv); CHKERRQ(ierr);
    n->nx = nv.x;
    n->ny = nv.y;
    n->nz = nv.z;

    // Tangential Vectors
    LocalCoor3DTangential( n );
    
    // Principal Curvatures
    j.x = GridBilinear( phi, GridFunction3D_DerivX, n->X );
    j.y = GridBilinear( phi, GridFunction3D_DerivY, n->X );
    j.z = GridBilinear( phi, GridFunction3D_DerivZ, n->X );
    IIMGlobalToLocal_1st( n, &j); // inverted rotation! global (x,y,z) to local (e,n,t) coor sys! 
    j.xx = GridBilinear( phi, GridFunction3D_DerivXX, n->X );
    j.xy = GridBilinear( phi, GridFunction3D_DerivXY, n->X );
    j.xz = GridBilinear( phi, GridFunction3D_DerivXZ, n->X );
    j.yy = GridBilinear( phi, GridFunction3D_DerivYY, n->X );
    j.yz = GridBilinear( phi, GridFunction3D_DerivYZ, n->X );
    j.zz = GridBilinear( phi, GridFunction3D_DerivZZ, n->X );
    IIMGlobalToLocal_2nd( n, &j); 
    n->k_nn = -j.nn / j.e;
    n->k_tt = -j.tt / j.e;
    n->k_nt = -j.nt / j.e;
    n->k  = GridBilinear( phi, GridFunction3D_Curv, n->X );
    //TODO: test curvatures values
    iim->F( n, iim->context );
  }
  //TODO: try to calculate # flops here
  ierr = PetscLogEventEnd(EVENT_IIMUpdateSurfaceQuantities,0,0,0,0); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
