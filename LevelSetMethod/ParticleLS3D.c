#include "LevelSetMethod.h"
#include "LSM_private.h"

#undef  __FUNCT__
#define __FUNCT__ "ParticleLS_AdvectParticles3D"
PetscErrorCode ParticleLS_AdvectParticles3D(ParticleLS pls, Coor dh, Grid velgrid, PetscReal dt)
{
  int i;
  Coor V1, V2, X;
  Particle p;
  PetscReal ****vel;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscLogEventBegin(EVENT_ParticleLS_AdvectParticles,0,0,0,0); CHKERRQ(ierr);
  ierr = GridGet(velgrid,&vel); CHKERRQ(ierr);
  for ( i = 0; i < ArrayLength(pls->particles); ++i) {
    ierr = ArrayGet(pls->particles,i,(void*)&p); CHKERRQ(ierr);
    /* RK2 Midpoint rule
     * V1 = U(Xp)
     * X  = Xp + dt/2 V1
     * V2 = U(X)
     * Xp = Xp + dt V2
     */
    ierr = InterpolateVelocity3D(0,vel,p->X,&V1); CHKERRQ(ierr);
    X.x = p->X.x + dt/2 * V1.x / dh.x;
    X.y = p->X.y + dt/2 * V1.y / dh.y;
    X.z = p->X.z + dt/2 * V1.z / dh.z;
    ierr = InterpolateVelocity3D(0,vel,X, &V2); CHKERRQ(ierr);
    p->X.x += dt * V2.x / dh.x;
    p->X.y += dt * V2.y / dh.y;
    p->X.z += dt * V2.z / dh.z;
  } // particles p
  ierr = PetscLogEventEnd(EVENT_ParticleLS_AdvectParticles,0,0,0,0); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef  __FUNCT__
#define __FUNCT__ "ParticleLS_ErrorCorrection3D"
PetscErrorCode ParticleLS_ErrorCorrection3D( ParticleLS pls, LevelSet ls )
{
  Particle p;
  PetscReal phi_p; // particle's phi value at grid points
  PetscReal dist; // particle's interpolated distance from interface
  iCoor A; // lower corner of the cell particle is located in
  iCoor a; // corner of cell in loop
  iCoor *b; // narrow band index
  int m,i,j,k;
  PetscReal ***phi_pos, ***phi_neg, ***phi;
  PetscReal sign;
  static Grid phi_pos_grid=0, phi_neg_grid=0;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscLogEventBegin(EVENT_ParticleLS_ErrorCorrection,0,0,0,0); CHKERRQ(ierr);
  if( phi_pos_grid == NULL ) {
    ierr = GridCreate(ls->phi->d,ls->phi->p,ls->phi->n,1,&phi_pos_grid); CHKERRQ(ierr);
    ierr = GridCreate(ls->phi->d,ls->phi->p,ls->phi->n,1,&phi_neg_grid); CHKERRQ(ierr);
  }
  ierr = GridCopy(ls->phi,phi_pos_grid); CHKERRQ(ierr); // phi+ = phi
  ierr = GridCopy(ls->phi,phi_neg_grid); CHKERRQ(ierr); // phi- = phi
  ierr = GridGet(phi_pos_grid,&phi_pos); CHKERRQ(ierr);
  ierr = GridGet(phi_neg_grid,&phi_neg); CHKERRQ(ierr);
  ierr = GridGet(ls->phi,&phi); CHKERRQ(ierr);
  for (m = 0; m < ArrayLength(pls->particles); ++m) {
    ierr = ArrayGet(pls->particles,m,&p); CHKERRQ(ierr);
    ierr = GridInterpolate( ls->phi, p->X, &dist ); CHKERRQ(ierr);
    sign = PetscSign(p->radius);
    // error identification
    // if a particle is classified as "escaped"
    if( sign * dist      < 0 &&        // particle on wrong side of interface
        PetscAbs( dist ) > PetscAbs( p->radius ) ) // by more than its radius
    {
      A.x = floor(p->X.x);
      A.y = floor(p->X.y);
      A.z = floor(p->X.y);

      // for each corner of the cell containing the particle
      for (k = 0; k < 2; ++k) {
        for (j = 0; j < 2; ++j) {
          for (i = 0; i < 2; ++i) {
            a.x = A.x + i;
            a.y = A.y + j;
            a.z = A.z + k;
            // quantify error: phi_p = Sp ( rp - | X - Xp | )
            phi_p = sign * ( sign*p->radius - sqrt( PetscSqr( p->X.x - a.x ) +
                                                    PetscSqr( p->X.y - a.y ) +
                                                    PetscSqr( p->X.z - a.z ) ) );
            if( sign > 0 ) {
              phi_pos[a.z][a.y][a.x] = PetscMax( phi_p, phi_pos[a.z][a.y][a.x] );
            } else {
              phi_neg[a.z][a.y][a.x] = PetscMin( phi_p, phi_neg[a.z][a.y][a.x] );
            } // phi+ or phi-
          } // i
        } // j
      } // k
    } // if escaped
  } // for each particle

  /* Merge phi+ and phi- into a single level set phi
   * phi = phi+ if |phi+| <= |phi-|
   *       phi- if |phi+| >  |phi-|
   */
  for( m = 0; m < ArrayLength(ls->band); m++ )
  {
    ierr = ArrayGet(ls->band,m,(void*)&b); CHKERRQ(ierr);
    if( PetscAbs(phi_pos[b->z][b->y][b->x]) <=
        PetscAbs(phi_neg[b->z][b->y][b->x]) ) {
      phi[b->z][b->y][b->x] = phi_pos[b->z][b->y][b->x];
    } else {
      phi[b->z][b->y][b->x] = phi_neg[b->z][b->y][b->x];
    }
  }
  ierr = PetscLogEventEnd(EVENT_ParticleLS_ErrorCorrection,0,0,0,0); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
