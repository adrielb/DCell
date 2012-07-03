#include "FluidField.h"
#include "FluidField_private.h"

PetscErrorCode FluidField_DiscreteCompatibilityCondition( FluidField f );

#undef __FUNCT__
#define __FUNCT__ "FluidFieldSolve"
PetscErrorCode FluidFieldSolve( FluidField f )
{
  PetscLogDouble t1,t2;
  PetscErrorCode ierr;
  PetscFunctionBegin;

  // AdvectSL each Eij
//  AdvectSL()

  // Rotate strain tensor by curl of velocity field
//  ???

  // Integrate strain rate
//  FluidFieldIntegrateStrainRate( DA daV, Vec vecV, DA daE, Vec vecE, PetscReal dh, PetscReal dt );

  // Compute the divergence of strain E->rhs
//  FluidFieldElasticDivergence()

  // IIM update for each level set

  // Enforce Dirichlet BC in case the IIM routines add correction terms to the boundaries
  ierr = FluidField_EnforceNoSlipBC(f); CHKERRQ(ierr);

  // Discrete compatibility condition
  ierr = FluidField_DiscreteCompatibilityCondition( f ); CHKERRQ(ierr);

  // Solve the system L(u) - px = del(sigma) + IIM
  ierr = PetscGetTime(&t1); CHKERRQ(ierr);
  ierr = KSPSolve(f->ksp,f->rhs,f->vel); CHKERRQ(ierr);
  ierr = PetscGetTime(&t2); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"KSP Solve: %f sec\n",t2-t1); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "FluidField_DiscreteCompatibilityCondition"
PetscErrorCode FluidField_DiscreteCompatibilityCondition( FluidField f )
{
  int i, lo, hi, bs;
  int count = 0;
  PetscReal avg, sum=0, eps = PETSC_MACHINE_EPSILON;
  PetscReal *rhs;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscLogEventBegin(EVENT_FluidField_DiscreteCompatibilityCondition,0,0,0,0); CHKERRQ(ierr);
  ierr = VecGetOwnershipRange(f->rhs, &lo, &hi); CHKERRQ(ierr);
  ierr = VecGetBlockSize(f->rhs, &bs); CHKERRQ(ierr);
  ierr = VecGetArray(f->rhs,&rhs); CHKERRQ(ierr);

  for ( i = lo+CELL_CENTER; i < hi; i+=bs ) {
    if( PetscAbs( rhs[i] ) > eps ) {
      sum += rhs[i];
      count++;
    }
  }

  ierr = MPI_Allreduce(&sum,&sum,1,MPI_DOUBLE,MPI_SUM,f->comm); CHKERRQ(ierr);
  ierr = MPI_Allreduce(&count,&count,1,MPI_INT,MPI_SUM,f->comm); CHKERRQ(ierr);
  avg = sum / count;
  ierr = PetscInfo1( 0, "avg C{div.u} = %e\n", avg); CHKERRQ(ierr);

  for ( i = lo+CELL_CENTER; i < hi; i+=bs ) {
    if( PetscAbs( rhs[i] ) > eps )
      rhs[i] -= avg;
  }

  ierr = VecRestoreArray(f->rhs,&rhs); CHKERRQ(ierr);
  ierr = PetscLogEventEnd(EVENT_FluidField_DiscreteCompatibilityCondition,0,0,0,0); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "FluidFieldMaxVelocityMag"
PetscErrorCode FluidFieldMaxVelocityMag( FluidField f, PetscReal *maxVel )
{
  iCoor e;
  Coor X,V;
  PetscReal *vel;
  PetscReal mag;
  int xs, ys, zs;
  int xm, ym, zm;
  DMDALocalInfo info;
  Vec locvel;
  DM dm = f->daV;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscLogEventBegin(EVENT_FluidFieldMaxVelocityMag,0,0,0,0); CHKERRQ(ierr);
  ierr = DMDAGetLocalInfo(dm,&info); CHKERRQ(ierr);
  ierr = DMDAGetCorners(dm,&xs,&ys,&zs,&xm,&ym,&zm); CHKERRQ(ierr);
  ierr = DMGetLocalVector(dm, &locvel); CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(dm, f->vel, INSERT_VALUES, locvel); CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(dm, f->vel, INSERT_VALUES, locvel); CHKERRQ(ierr);
  ierr = DMDAVecGetArrayDOF(dm,locvel,&vel); CHKERRQ(ierr);

  *maxVel = 0;
  e.x = xs+xm;
  e.y = ys+ym;
  e.z = zs+zm;
  if( e.x == info.mx ) e.x--;
  if( e.y == info.my ) e.y--;
  if( e.z == info.mz ) e.z--;

  if( f->is3D ) {
    PetscReal ****vel3D = (PetscReal****)vel;
    for ( X.z = zs; X.z < e.z; ++X.z) {
      for ( X.y = ys; X.y < e.y; ++X.y) {
        for ( X.x = xs; X.x < e.x; ++X.x) {
          ierr = InterpolateVelocity3D( U_FACE, vel3D, X, &V ); CHKERRQ(ierr);
          mag = PetscSqrtScalar( V.x*V.x + V.y*V.y + V.z*V.z );
          *maxVel = mag > *maxVel ? mag : *maxVel;
        }
      }
    }
  } else {
    PetscReal ***vel2D = (PetscReal***)vel;
    X.z = 0; // makes the compiler happy
    for ( X.y = ys; X.y < e.y; ++X.y) {
      for ( X.x = xs; X.x < e.x; ++X.x) {
        ierr = InterpolateVelocity2D( U_FACE, vel2D, X, &V ); CHKERRQ(ierr);
        mag = PetscSqrtScalar( V.x*V.x + V.y*V.y );
        *maxVel = mag > *maxVel ? mag : *maxVel;
      }
    }
  }

  ierr = MPI_Allreduce(maxVel,maxVel,1,MPI_DOUBLE,MPI_MAX,f->comm); CHKERRQ(ierr);

  ierr = DMDAVecRestoreArrayDOF(dm,f->vel,&vel); CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(dm,&locvel); CHKERRQ(ierr);
  ierr = PetscInfo1( 0, "max vel = %e\n", *maxVel); CHKERRQ(ierr);
  ierr = PetscLogEventEnd(EVENT_FluidFieldMaxVelocityMag,0,0,0,0); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
