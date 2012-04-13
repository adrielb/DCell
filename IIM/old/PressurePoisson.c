typedef struct _PressurePoisson
{
  IIM iim;
  int d1=63, d2=d1;
  Grid2D rhsC, p, px, py, u, v;
  LevelSet2D ls, lstemp;
  Mat m;
  KSP ksp;
} *PressurePoisson;

#undef __FUNCT__
#define __FUNCT__ "PressurePoissonCreate"
PetscInt EVENT_PressurePoissonCreate;
PetscErrorCode PressurePoissonCreate(  )
{
  PetscErrorCode ierr;
  
  
  PetscFunctionBegin;
  PetscLogEventBegin(EVENT_PressurePoissonCreate,0,0,0,0);
//  PetscLogEventRegister(&EVENT_PressurePoissonCreate,"PressurePoissonCreate", 0);
  ierr = KSPCreate(PETSC_COMM_SELF, &ksp); CHKERRQ(ierr);
  KSPGetPC(ksp, &pc);
//  KSPSetType(ksp, KSPCG);
//  PCSetType(pc, PCICC);
//  PCFactorSetLevels(pc, 0);
//  KSPSetInitialGuessNonzero(ksp, PETSC_TRUE);
  KSPSetType(ksp, KSPPREONLY);
  PCSetType(pc, PCCHOLESKY);
  PCFactorSetMatOrderingType(pc, MATORDERING_ND);
  KSPSetFromOptions(ksp);
  KSPSetOperators(ksp, m, m, SAME_PRECONDITIONER);
  PCSetUp(pc);
  
  ierr = IIMCreate(&iim, 12); CHKERRQ(ierr);
//  IIMSetForceComponents(iim, ForceComponentNormalSimple, ForceComponentTangentialSimple);
  CreateGrid2D(d1, d2, &rhsC);
  CreateGrid2D(d1, d2, &p);
  CreateGrid2D(d1, d2, &px);
  CreateGrid2D(d1, d2, &py);
  CreateGrid2D(d1, d2, &u);
  CreateGrid2D(d1, d2, &v);
  ierr = CreateLevelSet2D(d1,d2,&ls); CHKERRQ(ierr);  
  CreateLevelSet2D(d1, d2, &lstemp);
  LevelSetInitializeToStar(ls,PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE);
//  LevelSetInitializeToCircle(ls, PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE);
  
  PetscLogEventEnd(EVENT_PressurePoissonCreate,0,0,0,0);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PressurePoissonDetroy"
PetscInt EVENT_PressurePoissonDetroy;
PetscErrorCode PressurePoissonDetroy(  )
{
  PetscErrorCode ierr;
  
  
  PetscFunctionBegin;
  PetscLogEventBegin(EVENT_PressurePoissonDetroy,0,0,0,0);
//  PetscLogEventRegister(&EVENT_PressurePoissonDetroy,"PressurePoissonDetroy", 0);
  
  ierr = VecDestroy(magVel); CHKERRQ(ierr);
  ierr = VecDestroy(u2); CHKERRQ(ierr);
  ierr = VecDestroy(v2); CHKERRQ(ierr);
  ierr = DestroyLevelSet2D(ls); CHKERRQ(ierr);
  ierr = DestroyGrid2D(rhsC); CHKERRQ(ierr);
  ierr = DestroyGrid2D(p); CHKERRQ(ierr);
  ierr = DestroyGrid2D(px); CHKERRQ(ierr);
  ierr = DestroyGrid2D(py); CHKERRQ(ierr);
  ierr = DestroyGrid2D(u); CHKERRQ(ierr);
  ierr = DestroyGrid2D(v); CHKERRQ(ierr);
  ierr = IIMDestroy(iim); CHKERRQ(ierr);
  
  PetscLogEventEnd(EVENT_PressurePoissonDetroy,0,0,0,0);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PressurePoissonUpdate"
PetscInt EVENT_PressurePoissonUpdate;
PetscErrorCode PressurePoissonUpdate(  )
{
  PetscErrorCode ierr;
  
  
  PetscFunctionBegin;
  PetscLogEventBegin(EVENT_PressurePoissonUpdate,0,0,0,0);
//  PetscLogEventRegister(&EVENT_PressurePoissonUpdate,"PressurePoissonUpdate", 0);
  VecSet(rhsC->v, 0.); // TODO: this is important when reusing vecs, need to reset back to zero after time step
  ierr = IIMComputeCorrection2D(iim,JumpConditionPressure,ls,rhsC); CHKERRQ(ierr);
  IIMDiscreteCompatabilityCondition(ls, rhsC);
  ierr = KSPSolve(ksp, rhsC->v, p->v); CHKERRQ(ierr);
  IrregularNodeListWrite(ls,t);
  
  ierr = IIMPressureGradient(ls,p,px,py); CHKERRQ(ierr);
  
  ierr = IIMComputeCorrection2D(iim,JumpConditionXVelocity,ls,px); CHKERRQ(ierr);
  px->v2[d1/2][d2/2] = 0.;
  ierr = KSPSolve(ksp, px->v, u->v); CHKERRQ(ierr);
  
  ierr = IIMComputeCorrection2D(iim,JumpConditionYVelocity,ls,py); CHKERRQ(ierr);
  py->v2[d1/2][d2/2] = 0.;
  ierr = KSPSolve(ksp, py->v, v->v); CHKERRQ(ierr);
   
  VecPointwiseMult(u2, u->v, u->v);
  VecPointwiseMult(v2, u->v, u->v);
  VecWAXPY(magVel, 1, u2, v2);
  VecSqrt( magVel );

  VecNorm(magVel, NORM_INFINITY, &maxVel);
  dt = PetscMin(1/(2*maxVel),.1);
  LevelSetAdvect2D(dt, u, v, ls, lstemp);
  printf("\tmaxVel: %f",maxVel);
  printf("\tdt: %f\n",dt);

   
  
  PetscLogEventEnd(EVENT_PressurePoissonUpdate,0,0,0,0);
  PetscFunctionReturn(0);
}