#include "LevelSetMethod.h"

/* 
 * grad(F) . grad(phi) == 0
 * 
 */

PetscErrorCode LevelSetExtendVn_2D( LevelSet ls );
PetscErrorCode LevelSetExtendVn_3D( LevelSet ls );

#undef __FUNCT__
#define __FUNCT__ "LevelSetExtendVn"
PetscLogEvent EVENT_LevelSetExtendVn;
PetscErrorCode LevelSetExtendVn( LevelSet ls )
{
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  ierr = PetscLogEventBegin(EVENT_LevelSetExtendVn,0,0,0,0); CHKERRQ(ierr);
//  PetscLogEventRegister("LevelSetExtendVn", 0, &EVENT_LevelSetExtendVn);

  ierr = VecZeroEntries(ls->Vn->v); CHKERRQ(ierr);
  
  if( ls->phi->is2D )
  {
    ierr = LevelSetExtendVn_2D(ls); CHKERRQ(ierr);
  } else {
    ierr = LevelSetExtendVn_3D(ls); CHKERRQ(ierr);
  }
  
  PetscLogFlops( 0 );
  ierr = PetscLogEventEnd(EVENT_LevelSetExtendVn,0,0,0,0); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "LevelSetExtendVn_2D"
PetscErrorCode LevelSetExtendVn_2D( LevelSet ls )
{
  int i, j, b;
  PetscReal p, s, pxl, pxr, pyl, pyr, pxa, pxb, pya, pyb, sxa, sxb, sya, syb;
  Coor va, vb;
  PetscReal **phi = ls->phi->v2,
            **Vn  = ls->Vn->v2;
  iCoor *band;
  IrregularNode *n;
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  for( i = 0; i < ArrayLength(ls->irregularNodes); ++i) //Set the boundary conditions in Vn
  {
    ierr = ArrayGet(ls->irregularNodes,i,(void*)&n); CHKERRQ(ierr);
    Vn[n->y][n->x] = n->Vn;
  }
  
  for( b = 0; b < ArrayLength(ls->band); ++b)
  {
    ierr = ArrayGet(ls->band,b,(void*)&band); CHKERRQ(ierr);
    i = band->x;
    j = band->y;
    
    /*   yb
     *    |
     * xa-+-xb    xa--xl-+-xr--xb
     *    |
     *   ya
     */
//    ierr = VecWrite(ls->Vn->v,"iVn",b); CHKERRQ(ierr);
    if( PetscAbs(Vn[j][i]) > .0001 ) continue; //TODO: this test is not robust, need better way to prevent modifying boundary conditions

    p   = phi[j][i];
    s = PetscSign( p ); //Sign keeps track of upwinding 
    p = s*p;            //Invert the local stencil 
    pxl = s*phi[j][i-1];
    pxr = s*phi[j][i+1];
    pyl = s*phi[j-1][i];
    pyr = s*phi[j+1][i];
    
    pxa = p - pxl;  // D-phi 
    pxb = pxr - p;  // D+phi
    pya = p - pyl;  // D-phi
    pyb = pyr - p;  // D+phi
    
    sxa = pxa < 0 ? 0 :  1; // max( D-phi, 0 ) 
    sxb = pxb > 0 ? 0 : -1; // min( D+phi, 0 )
    sya = pya < 0 ? 0 :  1; // max( D-phi, 0 ) 
    syb = pyb > 0 ? 0 : -1; // min( D+phi, 0 )
    
    va.x = Vn[j][i-1];
    vb.x = Vn[j][i+1];
    va.y = Vn[j-1][i];
    vb.y = Vn[j+1][i];
    
    Vn[j][i] = ( va.x * pxa * sxa + vb.x * pxb * sxb + va.y * pya * sya + vb.y * pyb * syb ) /
                      ( pxa * sxa +        pxb * sxb +        pya * sya +        pyb * syb );
    
    if( Vn[j][i] != Vn[j][i] )
    {
      Vn[j][i] = 10;
      //TODO: make a better error handler for this NaN extension
      printf("Vn error: %f \n", ( pxa * sxa +        pxb * sxb +        pya * sya +        pyb * syb ) );
      exit(0);
    }
  }
  PetscLogFlops( 0 );
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "LevelSetExtendVn_3D"
PetscLogEvent EVENT_LevelSetExtendVn_3D;
PetscErrorCode LevelSetExtendVn_3D( LevelSet ls )
{
  int i, j, k, b;
  PetscReal p, s, pxl, pxr, pyl, pyr, pzl, pzr, 
                  pxa, pxb, pya, pyb, pza, pzb,
                  sxa, sxb, sya, syb, sza, szb;
  iCoor *band;
  Coor va, vb;
  PetscReal ***phi = ls->phi->v3;
  PetscReal ***Vn = ls->Vn->v3;
  IrregularNode *n;
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  ierr = PetscLogEventBegin(EVENT_LevelSetExtendVn_3D,0,0,0,0); CHKERRQ(ierr);
//  PetscLogEventRegister("LevelSetExtendVn_3D", 0, &EVENT_LevelSetExtendVn_3D);
  for( i = 0; i < ArrayLength(ls->irregularNodes); ++i) //Set the boundary conditions in Vn
  {
    ierr = ArrayGet(ls->irregularNodes,i,(void*)&n); CHKERRQ(ierr);
    Vn[n->z][n->y][n->x] = n->Vn;
  }
  
  for( b = 0; b < ArrayLength(ls->band); ++b)
  {
    ierr = ArrayGet(ls->band,b,(void*)&band); CHKERRQ(ierr);
    i = band->x;
    j = band->y;
    k = band->z;
    
//    if( PetscAbs(Vn[k][j][i]) > .0001 ) continue; //TODO: this test for BC node is not robust, need better way to skip boundary conditions from narrow band
    if( Vn[k][j][i] != 0 ) continue; 
    
    p   = phi[k][j][i];
    s = PetscSign( p ); //Sign keeps track of upwinding 
    p = s*p;            //Invert the local stencil
    
    pxl = s*phi[k][j][i-1];
    pxr = s*phi[k][j][i+1];
    pyl = s*phi[k][j-1][i];   // 6 flops?
    pyr = s*phi[k][j+1][i];
    pzl = s*phi[k-1][j][i];
    pzr = s*phi[k+1][j][i];
    
    pxa = p - pxl;  // D-phi 
    pxb = pxr - p;  // D+phi
    pya = p - pyl;  // D-phi  // 6 flops?
    pyb = pyr - p;  // D+phi
    pza = p - pzl;  // D-phi
    pzb = pzr - p;  // D+phi
    
    sxa = pxa < 0 ? 0 :  1; // max( D-phi, 0 ) 
    sxb = pxb > 0 ? 0 : -1; // min( D+phi, 0 )
    sya = pya < 0 ? 0 :  1; // max( D-phi, 0 ) // 6 flops? 
    syb = pyb > 0 ? 0 : -1; // min( D+phi, 0 )
    sza = pza < 0 ? 0 :  1; // max( D-phi, 0 ) 
    szb = pzb > 0 ? 0 : -1; // min( D+phi, 0 )
    
    va.x = Vn[k][j][i-1];
    vb.x = Vn[k][j][i+1];
    va.y = Vn[k][j-1][i];
    vb.y = Vn[k][j+1][i];
    va.z = Vn[k-1][j][i];
    vb.z = Vn[k+1][j][i];
    
    // 28 flops
    Vn[k][j][i] = ( va.x * pxa * sxa + vb.x * pxb * sxb + va.y * pya * sya + vb.y * pyb * syb + va.z * pza * sza + vb.z * pzb * szb ) /
                         ( pxa * sxa +        pxb * sxb +        pya * sya +        pyb * syb +        pza * sza +        pzb * szb );
  }
  
  PetscLogFlops( 0 );
  ierr = PetscLogEventEnd(EVENT_LevelSetExtendVn_3D,0,0,0,0); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
