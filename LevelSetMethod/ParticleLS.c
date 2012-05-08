#include "LevelSetMethod.h"
#include "LSM_private.h"

#undef __FUNCT__
#define __FUNCT__ "LevelSetInitializeParticles"
PetscErrorCode LevelSetInitializeParticles( LevelSet ls )
{
  int m,b;
  PetscReal dist;
  PetscReal phi_goal;
  PetscReal **phi2D, ***phi3D;
  ParticleLS pls;
  Particle p;
  iCoor *band;
  PetscReal sign;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscNew(struct _ParticleLS, &ls->pls); CHKERRQ(ierr);
  pls = ls->pls;
  pls->R_MIN = 0.1;
  pls->R_MAX = 0.5;
  pls->D_INIT= 3.0;
  pls->S_INIT = ls->phi->is2D ? 16 : 32;
  pls->L_MAX = 30;
  pls->G_TOL = 1e-3;

// use command line parameters for particle properties
  ierr = PetscOptionsGetReal(0,"-pls_rmin",&pls->R_MIN,0); CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(0,"-pls_rmax",&pls->R_MAX,0); CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(0,"-pls_dinit",&pls->D_INIT,0); CHKERRQ(ierr);
  ierr = PetscOptionsGetInt( 0,"-pls_sinit",&pls->S_INIT,0); CHKERRQ(ierr);

  ierr = ArrayCreate( "particles", sizeof(struct _Particle), &pls->particles); CHKERRQ(ierr);
  ierr = PetscRandomCreate(PETSC_COMM_WORLD,&pls->rnd); CHKERRQ(ierr);
  ierr = PetscRandomSetType(pls->rnd,PETSCRAND48); CHKERRQ(ierr);

  ierr = GridGet(ls->phi,&phi2D); CHKERRQ(ierr);
  ierr = GridGet(ls->phi,&phi3D); CHKERRQ(ierr);
  for ( b = 0; b < ArrayLength(ls->band); ++b) {
    ierr = ArrayGet(ls->band,b,(void*)&band); CHKERRQ(ierr);
    dist = ls->phi->is2D ? PetscAbs(phi2D[band->y][band->x]) :
                           PetscAbs(phi3D[band->z][band->y][band->x]);
    // if -D < phi(Xp) < D
    if( -pls->D_INIT < dist && dist < pls->D_INIT )
    {
      // number of particles seeded in each cell
      for ( m = 0; m < pls->S_INIT; ++m)
      {
        ierr = ArrayAppend(pls->particles,(void*)&p); CHKERRQ(ierr);

        PetscRandomGetValue(pls->rnd,&p->X.x);
        PetscRandomGetValue(pls->rnd,&p->X.y);
        PetscRandomGetValue(pls->rnd,&p->X.z);

        p->X.x += band->x - 0.5;
        p->X.y += band->y - 0.5;
        p->X.z += band->z - 0.5;

        PetscRandomGetValue(pls->rnd,&sign);
        sign = sign > 0.5 ? 1 : -1;

        PetscRandomGetValue(pls->rnd,&phi_goal);
        phi_goal = (pls->D_INIT - pls->R_MIN) * phi_goal + pls->R_MIN;
        phi_goal*= sign;

        p->radius = sign;
        ierr = ParticleLS_AttractParticleToPhiGoal( pls, ls, p, phi_goal ); CHKERRQ(ierr);
      } // m particles per cell
    } // if cell within +/-D of zero level set
  } // b in narrow band

  ierr = ParticleLS_AdjustRadii( pls, ls ); CHKERRQ(ierr);

  if( ls->phi->is2D ) {
    pls->AdvectParticles = ParticleLS_AdvectParticles2D;
    pls->ErrorCorrection = ParticleLS_ErrorCorrection2D;
  } else {
    pls->AdvectParticles = ParticleLS_AdvectParticles3D;
    pls->ErrorCorrection = ParticleLS_ErrorCorrection3D;
  }

  ls->Advect = LevelSetAdvectPLS;
  ierr = PetscInfo(0,"Particles initialized\n"); CHKERRQ(ierr);
  PetscFunctionReturn(0);
} // ParticleLSCreate()

#undef __FUNCT__
#define __FUNCT__ "ParticleLSDestroy"
PetscErrorCode ParticleLSDestroy( ParticleLS pls )
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = ArrayDestroy(pls->particles); CHKERRQ(ierr);
  if( pls->particles0 != NULL ) {
    ierr = ArrayDestroy(pls->particles0); CHKERRQ(ierr);
  }
  ierr = PetscRandomDestroy(&pls->rnd); CHKERRQ(ierr);
  ierr = PetscFree( pls ); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef  __FUNCT__
#define __FUNCT__ "ParticleLS_AttractParticleToPhiGoal"
PetscErrorCode ParticleLS_AttractParticleToPhiGoal( ParticleLS pls, LevelSet ls, Particle p, PetscReal phi_goal )
{
  Coor X = {0,0,0}; // New particle position
  Coor N; // normal direction at Xp
  PetscReal sign; // Sign( phi_goal - phi_Xp )
  PetscReal lambda = 1; // step size
  PetscReal phi_Xp, phi_X;
  Coor dh = ls->phi->d;
  PetscReal **phigrid;
  PetscErrorCode ierr;

  ierr = GridGet(ls->phi,&phigrid); CHKERRQ(ierr);

  int l;
  phi_Xp = Bilinear2D( GridFunction2D_Identity, phigrid, dh, p->X.x, p->X.y );
  for( l = 0; l < pls->L_MAX; ++l) {
    // Normal direction [N.x, N.y]
    N.x = Bilinear2D( GridFunction2D_DerivX, phigrid, dh, p->X.x, p->X.y );
    N.y = Bilinear2D( GridFunction2D_DerivY, phigrid, dh, p->X.x, p->X.y );
    PetscReal mag = sqrt( N.x*N.x + N.y*N.y );
    N.x = N.x / mag;
    N.y = N.y / mag;

    sign = PetscSign(phi_goal - phi_Xp);
    X.x = p->X.x + sign * lambda * N.x ;
    X.y = p->X.y + sign * lambda * N.y;
    phi_X = Bilinear2D( GridFunction2D_Identity, phigrid, dh, X.x, X.y );

    // step size failed, reduce it
    if( PetscAbs(phi_X - phi_goal) > PetscAbs(phi_Xp - phi_goal) ) {
      lambda = lambda / 2.;
    } else {
      p->X = X;
      phi_Xp = phi_X;
    }
    // if particle has reached phi_goal, break.
    if( PetscAbs(phi_X - phi_goal) < pls->G_TOL ) {
      break;
    }
  } // for l
  // if max iterations reached, particle may not have settled in appropriate band
  if( l == pls->L_MAX ) {
    phi_X = Bilinear2D( GridFunction2D_Identity, phigrid, dh, p->X.x, p->X.y);
    p->radius = PetscSign(phi_X);
  }
  return 0;
}

#undef  __FUNCT__
#define __FUNCT__ "ParticleLS_AdjustRadii"
PetscErrorCode ParticleLS_AdjustRadii( ParticleLS pls, LevelSet ls )
{
  int m;
  Particle p;
  PetscReal sphi;
  PetscReal dist;
  PetscReal sign;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscLogEventBegin(EVENT_ParticleLS_AdjustRadii,0,0,0,0); CHKERRQ(ierr);
  for ( m = 0; m < ArrayLength(pls->particles); ++m) {
    ierr = ArrayGet(pls->particles,m,&p); CHKERRQ(ierr);
    // assign particle radius
    ierr = GridInterpolate( ls->phi, p->X, &dist); CHKERRQ(ierr);
    sign = PetscSign(p->radius);
    sphi = sign * dist;
    if( pls->R_MIN < sphi && sphi < pls->R_MAX ) {
      p->radius = sign*sphi;
    } else if( sphi > pls->R_MAX ) {
      p->radius = sign * pls->R_MAX;
    } else if( sphi < pls->R_MIN ) {
      p->radius = sign * pls->R_MIN;
    } // radius
  } // particles
  ierr = PetscLogEventEnd(EVENT_ParticleLS_AdjustRadii,0,0,0,0); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


#undef  __FUNCT__
#define __FUNCT__ "ParticleLSWriteParticles"
PetscErrorCode ParticleLSWriteParticles(ParticleLS pls, int t)
{
  PetscErrorCode ierr;
  PetscFunctionBegin;
  if( pls ) {
    ierr = ArrayWrite(pls->particles, t); CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

#undef  __FUNCT__
#define __FUNCT__ "ParticleLS_ReseedParticles"
PetscErrorCode ParticleLS_ReseedParticles( ParticleLS pls, LevelSet ls )
{
  int i, b;
  iCoor *band, cell;
  PetscReal dist, sign;
  PetscReal **phi;
//  PetscReal phimin, phimax, phi_goal;
  PetscReal phi_goal;
  Particle p;
  int *count;
  const int targetDensity = pls->S_INIT;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  static Array particleCounts, postCounts;
  if( particleCounts == NULL ) {
    ierr = ArrayCreate( "particleCounts", sizeof(int), &particleCounts); CHKERRQ(ierr);
    ierr = ArrayCreate( "postCounts", sizeof(int), &postCounts); CHKERRQ(ierr);
  }
  ierr = PetscLogEventBegin(EVENT_ParticleLS_ReseedParticles,0,0,0,0); CHKERRQ(ierr);
  ierr = GridGet(ls->phi,&phi); CHKERRQ(ierr);
  /* Delete particles that are
   *   1) too far: p->sign*phi(p->X) > D_INIT
   *   2) TODO: too dense: delete particles farthest from the interface
   */
  for ( i = 0; i < ArrayLength(pls->particles); ++i) {
    ierr = ArrayGet(pls->particles,i,&p); CHKERRQ(ierr);

    ierr = GridInterpolate( ls->phi, p->X, &dist ); CHKERRQ(ierr);
    if( PetscAbs(dist) > pls->D_INIT ) {
      // TODO: this deletes both escaped and non-escaped particles. Is it possible for an escaped particle to be so far from an interface?
      // If so, it may not be in the level set bounding box and thus not have a velocity field for interpolation
      ierr = ArrayDelete1(pls->particles, i); CHKERRQ(ierr);
      i--;
      continue;
    }
  }

  // Count number of particles in each cell.
  ierr = ArraySetCoor(particleCounts,ls->phi->p,ls->phi->n); CHKERRQ(ierr);
  ierr = ArrayZero(particleCounts); CHKERRQ(ierr);
  for ( i = 0; i < ArrayLength(pls->particles); ++i) {
    ierr = ArrayGet(pls->particles,i,&p); CHKERRQ(ierr);
    ierr = GridInterpolate( ls->phi, p->X, &sign); CHKERRQ(ierr);
    sign = PetscSign(sign);

    cell.x = (int)(p->X.x + 0.5);
    cell.y = (int)(p->X.y + 0.5);
    cell.z = (int)(p->X.z + 0.5);

    ArrayGetCoor(particleCounts, cell, &count);
    (*count)++;

    // if its a non-escaped particle
    if(p->radius * sign > 0) {
    }
  }
  static int filecount;

  // Add particles where density too low
  for ( b = 0; b < ArrayLength(ls->band); ++b) {
    ierr = ArrayGet(ls->band,b,(void*)&band); CHKERRQ(ierr);
    dist = PetscAbs(phi[band->y][band->x]);
    if( dist > pls->D_INIT ) continue;
    ierr = ArrayGetCoor(particleCounts, *band, &count); CHKERRQ(ierr);
    if(  (*count) > 1e3 ) {
      printf("*****count: %d\n", *count );
      printf("filecount: %d\n", filecount );
      printf("%d: { %d, %d, %d}\n", b, band->x, band->y, band->z );
      exit(0);
    }

    if( filecount == 153 ) {
      int *c;
      iCoor a = {47, 0, 0};
      ierr = ArrayGetCoor(particleCounts, a, &c); CHKERRQ(ierr);
      printf("@@@@c: %d \t %d: {%d, %d} %d \n", *c, b, band->x, band->y, *count );
      iCoor A = {47,-1, 0};
      ierr = ArrayGetCoor(particleCounts, A, &c); CHKERRQ(ierr);
      printf("$$$$c: %d \t %d: {%d, %d} %d \n", *c, b, band->x, band->y, *count );
    }
    for (; *count < targetDensity; (*count)++) {
      ierr = ArrayAppend(pls->particles,(void*)&p); CHKERRQ(ierr);

      PetscRandomGetValue(pls->rnd,&p->X.x);
      PetscRandomGetValue(pls->rnd,&p->X.y);
      PetscRandomGetValue(pls->rnd,&p->X.z);

      p->X.x += band->x - 0.5;
      p->X.y += band->y - 0.5;
      p->X.z += band->z - 0.5;

      ierr = GridInterpolate( ls->phi, p->X, &phi_goal); CHKERRQ(ierr);
      sign = PetscSign(phi_goal);
      phi_goal = PetscAbs(phi_goal);
      phi_goal = PetscMax(pls->R_MIN,phi_goal);
      phi_goal = PetscMin(pls->D_INIT,phi_goal);
      phi_goal = phi_goal * sign;
/*
      PetscRandomGetValue(pls->rnd,&phi_goal);
      phimax = PetscMin(pls->R_MAX,dist+0.71);
      phimin = PetscMax(pls->R_MIN,dist-0.71);
      phi_goal = (phimax - phimin) * phi_goal + phimin;
      phi_goal*= sign;
*/
      p->radius = sign;
//      ierr = ParticleLS_AttractParticleToPhiGoal( pls, ls, p, phi_goal ); CHKERRQ(ierr);
    } // add particle to cell b
  } // for each b in band

    // Count number of particles in each cell.
  	ierr = ArraySetCoor(postCounts,ls->phi->p,ls->phi->n); CHKERRQ(ierr);
    ierr = ArrayZero(postCounts); CHKERRQ(ierr);
    for ( i = 0; i < ArrayLength(pls->particles); ++i) {
      ierr = ArrayGet(pls->particles,i,&p); CHKERRQ(ierr);
      cell.x = (int)(p->X.x + 0.5);
      cell.y = (int)(p->X.y + 0.5);
      cell.z = (int)(p->X.z + 0.5);
      ArrayGetCoor(postCounts, cell, &count);
      (*count)++;
    }

    ierr = ArrayWrite(particleCounts, filecount); CHKERRQ(ierr);
    ierr = ArrayWrite(postCounts, filecount); CHKERRQ(ierr);
    filecount++;

  ierr = PetscLogEventEnd(EVENT_ParticleLS_ReseedParticles,0,0,0,0); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef  __FUNCT__
#define __FUNCT__ "LevelSetAdvectPLS"
PetscErrorCode LevelSetAdvectPLS(LevelSet ls, Grid velgrid, PetscReal dt)
{
  ParticleLS pls = ls->pls;
  PetscErrorCode ierr;
  /* 1) Evolve particles by RK2
   * 2) Evolve LS by SL1
   * 3) Error correction
   * 4) Reinit LS w/FMM
   * 5) Error correction
   * 6) Adjust particle radii
   * 7) Reseed particles
   */

  PetscFunctionBegin;
  ierr = pls->AdvectParticles(pls,ls->phi->d,velgrid,dt); CHKERRQ(ierr);
  ierr = GridCopy(ls->phi,ls->phi0); CHKERRQ(ierr);
  ierr = LevelSetAdvectSL(ls,velgrid,dt); CHKERRQ(ierr);
  ierr = pls->ErrorCorrection(pls,ls); CHKERRQ(ierr);
  ierr = LevelSetUpdateIrregularNodeList( ls ); CHKERRQ(ierr);
  ierr = LevelSetCFLIncrement( ls, velgrid, dt ); CHKERRQ(ierr);
  if( ls->AdvectCount >= ls->AdvectThres || ls->CFLcount >= ls->CFLthres ) {
    ierr = PetscInfo4(0,"CFLcount: %f:%f  AdvectCount: %d:%d\n", ls->CFLcount, ls->CFLthres, ls->AdvectCount,ls->AdvectThres); CHKERRQ(ierr);
    ierr = LevelSetReinitialize(ls); CHKERRQ(ierr);
    ierr = pls->ErrorCorrection(pls,ls); CHKERRQ(ierr);
    ierr = ParticleLS_ReseedParticles(pls,ls); CHKERRQ(ierr);
    ierr = ParticleLS_AdjustRadii(pls,ls); CHKERRQ(ierr);
    ls->CFLcount = 0;
    ls->AdvectCount = 0;
  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "LevelSetAdvectPLSRK2HalfStep"
PetscErrorCode LevelSetAdvectPLSRK2HalfStep( LevelSet ls, Grid velgrid, PetscReal dt )
{
  ParticleLS pls = ls->pls;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  if( pls->particles0 == NULL ) {
    ierr = ArrayDuplicate(pls->particles,&pls->particles0); CHKERRQ(ierr);
  }
  ierr = ArrayCopy(pls->particles,pls->particles0); CHKERRQ(ierr);
  ierr = ParticleLS_AdvectParticles2D(pls,ls->phi->d,velgrid,dt); CHKERRQ(ierr);
  ierr = GridCopy(ls->phi,ls->phi0); CHKERRQ(ierr);
  ierr = LevelSetAdvectSL(ls, velgrid, dt ); CHKERRQ(ierr);
  ierr = pls->ErrorCorrection(pls,ls); CHKERRQ(ierr);
  ierr = LevelSetUpdateIrregularNodeList( ls ); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "LevelSetAdvectPLSRK2FullStep"
PetscErrorCode LevelSetAdvectPLSRK2FullStep( LevelSet ls, Grid velgrid, PetscReal dt )
{
  ParticleLS pls = ls->pls;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = ArrayCopy(pls->particles0,pls->particles); CHKERRQ(ierr);
  ierr = ParticleLS_AdvectParticles2D(pls,ls->phi->d,velgrid,dt); CHKERRQ(ierr);
  ierr = LevelSetAdvectSL(ls, velgrid, dt ); CHKERRQ(ierr);
  ierr = pls->ErrorCorrection(pls,ls); CHKERRQ(ierr);
  ierr = LevelSetUpdateIrregularNodeList( ls ); CHKERRQ(ierr);
  ierr = LevelSetCFLIncrement( ls, velgrid, dt ); CHKERRQ(ierr);
  if( ls->AdvectCount >= ls->AdvectThres || ls->CFLcount >= ls->CFLthres ) {
    ierr = PetscInfo4(0,"CFLcount: %f:%f  AdvectCount: %d:%d\n", ls->CFLcount, ls->CFLthres, ls->AdvectCount,ls->AdvectThres); CHKERRQ(ierr);
    ierr = LevelSetReinitialize(ls); CHKERRQ(ierr);
    ierr = pls->ErrorCorrection(pls,ls); CHKERRQ(ierr);
    ierr = ParticleLS_ReseedParticles(pls,ls); CHKERRQ(ierr);
    ierr = ParticleLS_AdjustRadii(pls,ls); CHKERRQ(ierr);
    ls->CFLcount = 0;
    ls->AdvectCount = 0;
  }
  PetscFunctionReturn(0);
}


/*
 * for multiple interacting level sets
 *   seed interior region
 *   rebuild phi- from LS particles
 *   rebuild phi+ from all other region particles
 *   use particles to correct after advection and after reinit
 *   apply projection method for each correction step
 *   as a result, two level sets are updated locally for each interface
 * projection method
 *   subtract the avg of the two smallest phi_i from all the other phi_i
 */
