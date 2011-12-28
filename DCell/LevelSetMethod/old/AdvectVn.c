#include "LevelSetMethod.h"

/*
 * phi_t + Vn |del(phi)| == 0
 */
PetscErrorCode LevelSetAdvectVn_2D( LevelSet ls, PetscReal dt );
PetscErrorCode LevelSetAdvectVn_3D( LevelSet ls, PetscReal dt );

#undef __FUNCT__
#define __FUNCT__ "LevelSetAdvectVn"
PetscInt EVENT_LevelSetAdvectVn;
PetscErrorCode LevelSetAdvectVn( LevelSet ls  )
{
  PetscReal maxVn=0, aVn, dt, dx;
  IrregularNode *n;
  Coor d = ls->Vn->d;
  int i;
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  ierr = PetscLogEventBegin(EVENT_LevelSetAdvectVn,0,0,0,0); CHKERRQ(ierr);
//  PetscLogEventRegister(&EVENT_LevelSetAdvectVn,"LevelSetAdvectVn", 0);
  ierr = LevelSetExtendVn( ls ); CHKERRQ(ierr);
  
  // Determine time step: dt |Vn| < dx
  for( i = 0; i < ArrayLength(ls->irregularNodes); i++)
  {
    ierr = ArrayGet(ls->irregularNodes, i,(void*)&n); CHKERRQ(ierr);
    aVn = PetscAbs(n->Vn);
    maxVn = maxVn < aVn ? aVn : maxVn;
  }
  PetscReal mindt = 1e-4;
  PetscReal maxCFL= 0.5;
  dx = PetscMin( d.x, d.y );
  dt = maxCFL * dx / maxVn;
  dt = PetscMin(mindt,dt);
  
  printf("\t dt: %1.2e \t maxDt: %1.2e \t CFL: %1.2e \n", dt, dx/maxVn, maxVn * dt / dx );

  if( ls->phi->is2D )
  {
    ierr = LevelSetAdvectVn_2D( ls, dt); CHKERRQ(ierr);
  } else {
    ierr = LevelSetAdvectVn_3D( ls, dt); CHKERRQ(ierr);
  }
  
  iCoor CELL_CENTER={0,0,0};
  ierr = IrregularNodeListUpdate(CELL_CENTER, ls); CHKERRQ(ierr);
  ierr = LevelSetReinitialize(ls); CHKERRQ(ierr);

  ierr = PetscLogEventEnd(EVENT_LevelSetAdvectVn,0,0,0,0); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "LevelSetAdvectVn_2D"
PetscErrorCode LevelSetAdvectVn_2D( LevelSet ls, PetscReal dt )
{
  int i,j;
  iCoor *band;
  PetscReal a,z;
  PetscReal magX, magY, mag;
  PetscReal **phi = ls->phi->v2,
            **Vn  = ls->Vn->v2,
            **fi  = ls->temp->v2;
  Coor d = ls->phi->d;
  int b;
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  ierr = VecCopy(ls->phi->v,ls->temp->v); CHKERRQ(ierr); //TODO: get rid of this veccopy
  // Determine |del(phi)| and advect phi
  for ( b = 0; b < ArrayLength(ls->band); ++b)
  {
    ierr = ArrayGet(ls->band,b,(void*)&band); CHKERRQ(ierr);
    i = band->x;
    j = band->y;
    
    /* |del(phi)| = ( max( s D-x, -s D+x, 0)^2 + max( s D-y, -s D+y, 0)^2 )^0.5
    s = PetscSign( Vn[j][i] );
    magX = PetscMax( s * (phi[j][i] - phi[j][i-1]) / d.x, 
                    -s * (phi[j][i+1] - phi[j][i]) / d.x );
    magX = PetscMax( magX, 0 );
    magY = PetscMax( s * (phi[j][i] - phi[j-1][i]) / d.y,
                    -s * (phi[j+1][i] - phi[j][i]) / d.y );
    magY = PetscMax( magY, 0 );
    */
    if( Vn[j][i] > 0 )
    {
      a = PetscMax( (phi[j][i] - phi[j][i-1]) / d.x, 0);
      z = PetscMin( (phi[j][i+1] - phi[j][i]) / d.x, 0);
      magX = PetscMax( a*a , z*z );

      a = PetscMax( (phi[j][i] - phi[j-1][i]) / d.y, 0);
      z = PetscMin( (phi[j+1][i] - phi[j][i]) / d.y, 0);
      magY = PetscMax( a*a , z*z );
    } else {
      a = PetscMin( (phi[j][i] - phi[j][i-1]) / d.x, 0);
      z = PetscMax( (phi[j][i+1] - phi[j][i]) / d.x, 0);
      magX = PetscMax( a*a , z*z );

      a = PetscMin( (phi[j][i] - phi[j-1][i]) / d.y, 0);
      z = PetscMax( (phi[j+1][i] - phi[j][i]) / d.y, 0);
      magY = PetscMax( a*a , z*z );
    }

    mag = sqrt( magX + magY );
    
    // phi_1 = phi_0 - Vn * |grad(phi_0)| * dt
    fi[j][i] = phi[j][i] - Vn[j][i] * mag * dt;
  }
  ierr = VecCopy(ls->temp->v,ls->phi->v); CHKERRQ(ierr);
  
  PetscLogFlops( 0 );
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "LevelSetAdvectVn_3D"
PetscErrorCode LevelSetAdvectVn_3D( LevelSet ls, PetscReal dt )
{
  int i,j,k;
  iCoor *band;
  PetscReal magX, magY, magZ, mag;
  PetscReal ***phi = ls->phi->v3,
            ***Vn  = ls->Vn->v3,
            ***fi  = ls->temp->v3;
  PetscReal s;
  Coor d = ls->phi->d;
  int b;
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  ierr = VecCopy(ls->phi->v,ls->temp->v); CHKERRQ(ierr); //TODO: this veccopy is needed to preserve topology in narrow band,
  
  // Determine |del(phi)| and advect phi
  for( b = 0; b < ArrayLength(ls->band); ++b)
  {
    ierr = ArrayGet(ls->band,b,(void*)&band); CHKERRQ(ierr);
    i = band->x;
    j = band->y;
    k = band->z;

    s = PetscSign( Vn[k][j][i] );
    
    // |del(phi)| = ( max( s D-x, -s D+x, 0)^2 + max( s D-y, -s D+y, 0)^2 )^0.5
    magX = PetscMax( s * (phi[k][j][i] - phi[k][j][i-1]) / d.x, 
                    -s * (phi[k][j][i+1] - phi[k][j][i]) / d.x );
    magX = PetscMax( magX, 0 );
    magY = PetscMax( s * (phi[k][j][i] - phi[k][j-1][i]) / d.y,
                    -s * (phi[k][j+1][i] - phi[k][j][i]) / d.y );
    magY = PetscMax( magY, 0 );
    magZ = PetscMax( s * (phi[k][j][i] - phi[k-1][j][i]) / d.z, 
                    -s * (phi[k+1][j][i] - phi[k][j][i]) / d.z );
    magZ = PetscMax( magZ, 0 );
    mag = sqrt( magX*magX + magY*magY + magZ*magZ );
    
    // phi_1 = phi_0 - Vn * |grad(phi)| * dt
    fi[k][j][i] = phi[k][j][i] - Vn[k][j][i] * mag * dt;
  }
  ierr = VecCopy(ls->temp->v,ls->phi->v); CHKERRQ(ierr);
  
  PetscLogFlops( 0 );
  PetscFunctionReturn(0);
}
