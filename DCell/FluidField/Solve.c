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

  // Discrete compatibility condition
  ierr = FluidField_DiscreteCompatibilityCondition( f ); CHKERRQ(ierr);

  // Enforce Dirichlet BC in case the IIM routines add correction terms to the boundaries
  ierr = FluidField_EnforceNoSlipBC(f); CHKERRQ(ierr);

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
  PetscReal avg, sum=0, eps = PETSC_MACHINE_EPSILON;
  PetscReal ***rhs;
  int i,j;
  int count = 0;
  int xs,ys,xm,ym;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscLogEventBegin(EVENT_FluidField_DiscreteCompatibilityCondition,0,0,0,0); CHKERRQ(ierr);
  ierr = DAGetCorners(f->daV,&xs,&ys,0,&xm,&ym,0); CHKERRQ(ierr);
  ierr = DAVecGetArrayDOF(f->daV,f->rhs,&rhs); CHKERRQ(ierr);
  for ( j = ys; j < ys+ym; ++j) {
    for ( i = xs; i < xs+xm; ++i) {
      if( PetscAbs(rhs[j][i][CELL_CENTER]) > eps )
      {
        sum += rhs[j][i][CELL_CENTER];
        count++;
      }
    }
  }
  ierr = MPI_Allreduce(&sum,&sum,1,MPI_DOUBLE,MPI_SUM,f->comm); CHKERRQ(ierr);
  ierr = MPI_Allreduce(&count,&count,1,MPI_INT,MPI_SUM,f->comm); CHKERRQ(ierr);
  avg = sum / count;
  for ( j = ys; j < ys+ym; ++j) {
    for ( i = xs; i < xs+xm; ++i) {
      if( PetscAbs(rhs[j][i][CELL_CENTER]) > eps )
      {
        rhs[j][i][CELL_CENTER] -= avg;
      }
    }
  }
  ierr = DAVecRestoreArrayDOF(f->daV,f->rhs,&rhs); CHKERRQ(ierr);
  ierr = PetscLogEventEnd(EVENT_FluidField_DiscreteCompatibilityCondition,0,0,0,0); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "FluidFieldMaxVelocityMag"
PetscErrorCode FluidFieldMaxVelocityMag( FluidField f, PetscReal *maxVel )
{
  iCoor e;
  Coor X,V;
  PetscReal ***vel;
  PetscReal mag;
  int xs,ys,xm,ym;
  DALocalInfo info;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscLogEventBegin(EVENT_FluidFieldMaxVelocityMag,0,0,0,0); CHKERRQ(ierr);
  ierr = DAGetLocalInfo(f->daV,&info); CHKERRQ(ierr);
  ierr = DAGetCorners(f->daV,&xs,&ys,0,&xm,&ym,0); CHKERRQ(ierr);
  ierr = DAVecGetArrayDOF(f->daV,f->vel,&vel); CHKERRQ(ierr);

  // hack: ignore upper boundary since vel interpolation needs [x, x+1]
  e.x = xs+xm-1;
  e.y = ys+ym-1;
  //TODO: Warning: doesnt count velocity field between proc boundaries, need to do
  //      global to local sync to include ghost nodes
  *maxVel = 0;
  for ( X.y = ys; X.y < e.y; ++X.y) {
    for ( X.x = xs; X.x < e.x; ++X.x) {
      ierr = InterpolateVelocity2D( U_FACE, vel, X, &V ); CHKERRQ(ierr);
      mag = PetscSqrtScalar( V.x*V.x + V.y*V.y );
      *maxVel = mag > *maxVel ? mag : *maxVel;
    }
  }
  ierr = DAVecRestoreArrayDOF(f->daV,f->vel,&vel); CHKERRQ(ierr);
  ierr = MPI_Allreduce(maxVel,maxVel,1,MPI_DOUBLE,MPI_MAX,f->comm); CHKERRQ(ierr);
  ierr = PetscLogEventEnd(EVENT_FluidFieldMaxVelocityMag,0,0,0,0); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
