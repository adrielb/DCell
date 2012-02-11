#ifndef LSM_PRIVATE_H_
#define LSM_PRIVATE_H_

PetscLogEvent EVENT_LevelSetUpdateIrregularNodeList;
PetscLogEvent EVENT_LevelSetReinitialize;
PetscLogEvent EVENT_LevelSetWriteIrregularNodeList;
PetscLogEvent EVENT_LevelSetGetVelocity;
PetscLogEvent EVENT_LevelSetAdvectSL;
PetscLogEvent EVENT_ParticleLS_AdvectParticles;
PetscLogEvent EVENT_ParticleLS_ErrorCorrection;
PetscLogEvent EVENT_ParticleLS_AdjustRadii;
PetscLogEvent EVENT_ParticleLS_ReseedParticles;

PetscErrorCode LevelSetRegisterEvents( );
PetscErrorCode LevelSetUpdateIrregularNodeList_2D( LevelSet ls, Grid p );
PetscErrorCode LevelSetUpdateIrregularNodeList_3D( LevelSet ls, Grid p );
PetscErrorCode LevelSetGetVelocity(LevelSet ls, int ga, Grid velgrid);
PetscErrorCode LevelSetGatherVelocity(LevelSet ls, int ga, Grid velgrid);
PetscErrorCode LevelSetAdvectSL(LevelSet ls, Grid velgrid, PetscReal dt);
PetscErrorCode LevelSetCFLIncrement( LevelSet ls, Grid velgrid, PetscReal dt );
PetscErrorCode OrthogonalProjection2D( double phi3[3][3], double phi[5][5], Coor *op );
PetscErrorCode OrthogonalProjection2D_Quadratic( double phi[3][3],  Coor *op );
PetscErrorCode OrthogonalProjection2D_Linear( double phi[3][3],  Coor *op );
PetscErrorCode OrthogonalProjection3D( double phi[3][3][3], Coor *op );
PetscErrorCode OrthogonalProjection3D_1st( double phi[3][3][3], Coor *op );

typedef struct _FMMNode
{
  PetscReal phi;
  iCoor pos;
} *FMMNode;

PetscErrorCode FMM_InitializeEikonal( LevelSet ls, MemCache mc, Heap heap );
PetscErrorCode FMM_SolveEikonal( LevelSet ls, MemCache mc, Heap heap, int sign );
PetscErrorCode FMM_PushNeighbors2D( LevelSet ls, MemCache mc, Heap heap, int sign, iCoor pos );
PetscErrorCode FMM_PushNeighbors3D( LevelSet ls, MemCache mc, Heap heap, int sign, iCoor pos );




typedef struct _Particle {
  Coor X;
  PetscReal radius;
} *Particle;

struct _ParticleLS {
  /* PARAMETERS */
  PetscReal R_MIN;   // minimum particle radius
  PetscReal R_MAX;   // maximum particle radius
  PetscReal D_INIT;  // initial placement with +/-3 of zero LS
  PetscReal S_INIT;    // initial seeding density
  int L_MAX;    // max iterations for finding phi_goal
  PetscReal G_TOL;   // tolerance when phi == phi_goal

  /* ACTIONS */
  PetscErrorCode (*AdvectParticles)(ParticleLS pls, Coor dh, Grid velgrid, PetscReal dt);
  PetscErrorCode (*ErrorCorrection)(ParticleLS pls, LevelSet ls);

  /* OBJECTS */
  PetscRandom rnd;   // U[0,1]
  Array particles;
  Array particles0; // particle locations at beginning of step in RK2
};

PetscErrorCode ParticleLS_AttractParticleToPhiGoal( ParticleLS pls, LevelSet ls, Particle p, PetscReal phi_goal );
PetscErrorCode ParticleLS_AdjustRadii( ParticleLS pls, LevelSet ls );
PetscErrorCode ParticleLS_AdvectParticles2D(ParticleLS pls, Coor dh, Grid velgrid, PetscReal dt);
PetscErrorCode ParticleLS_AdvectParticles3D(ParticleLS pls, Coor dh, Grid velgrid, PetscReal dt);
PetscErrorCode ParticleLS_ErrorCorrection2D( ParticleLS pls, LevelSet ls );
PetscErrorCode ParticleLS_ErrorCorrection3D( ParticleLS pls, LevelSet ls );

#endif /* LSM_PRIVATE_H_ */
