#include "FluidField.h"

PetscReal OptimalSOR( DA da );

#undef __FUNCT__
#define __FUNCT__ "FluidFieldCreate"
PetscInt EVENT_FluidFieldCreate;
PetscErrorCode FluidFieldCreate(FluidField *fluid)
{
  PC pc;
  FluidField f;
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  PetscLogEventBegin(EVENT_FluidFieldCreate,0,0,0,0);
//  PetscLogEventRegister(&EVENT_FluidFieldCreate,"FluidFieldCreate", 0);
  ierr = PetscNew(struct _FluidField, &f); CHKERRQ(ierr);
  /*
   *  Default values
   */ 
  f->mu = 1;
  f->rho = 1;
  PetscReal dx = 1;
  int MX = 64; //default dimensions if command line opts not specified
  
  /*TODO: implement two different ways for scaling length
  f->d.x = f->l.x / (i.mx-2);
  f->d.y = f->l.y / (i.my-2);
  f->d.z = f->l.z / (i.mz-2);
  */
  
  int dummy;
  ierr = PetscOptionsGetReal(0,"-fluid_rho",&f->rho,0); CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(0,"-fluid_mu",&f->mu,0); CHKERRQ(ierr);
  ierr = PetscOptionsGetInt(0,"-da_grid_z",&dummy,&f->is3D); CHKERRQ(ierr);
  if( f->is3D )
  {
    //TODO: try using MPI_Dims_Create to provide values for lx,ly,lz
    ierr = DACreate3d(PETSC_COMM_WORLD, DA_NONPERIODIC, DA_STENCIL_STAR,
          -MX,-MX,-MX,PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE,1,1,0,0,0,&f->da); CHKERRQ(ierr);
    ierr = FluidFieldAssembleStaggeredGrid_3D( dx, f); CHKERRQ(ierr);
  } else {
    ierr = DACreate2d(PETSC_COMM_WORLD,DA_NONPERIODIC,DA_STENCIL_STAR,
          -MX,-MX,PETSC_DECIDE,PETSC_DECIDE,1,1,0,0, &f->da); CHKERRQ(ierr);
    ierr = FluidFieldAssembleStaggeredGrid_2D( dx, f); CHKERRQ(ierr);
  }

  ierr = DACreateGlobalVector(f->da,&f->p); CHKERRQ(ierr);
  ierr = VecDuplicate(f->p,&f->source); CHKERRQ(ierr);
  ierr = DACreateGlobalArray(f->da,&f->gaU,&f->u); CHKERRQ(ierr);
  ierr = DACreateGlobalArray(f->da,&f->gaV,&f->v); CHKERRQ(ierr);
//  ierr = VecDuplicate(f->p,&f->u); CHKERRQ(ierr);
//  ierr = VecDuplicate(f->p,&f->v); CHKERRQ(ierr);
  ierr = VecDuplicate(f->p,&f->px); CHKERRQ(ierr);
  ierr = VecDuplicate(f->p,&f->py); CHKERRQ(ierr);
  ierr = VecDuplicate(f->p,&f->div); CHKERRQ(ierr);

  if( f->is3D )
  {
    ierr = DACreateGlobalArray(f->da,&f->gaW,&f->w); CHKERRQ(ierr);
//    ierr = VecDuplicate(f->p,&f->w); CHKERRQ(ierr);
    ierr = VecDuplicate(f->p,&f->pz); CHKERRQ(ierr);
  }
  
  /* Velocity Solver Setup */
  ierr = KSPCreate(PETSC_COMM_WORLD,&f->kspU); CHKERRQ(ierr);
  ierr = KSPSetOperators(f->kspU,f->matU,f->matU,SAME_PRECONDITIONER); CHKERRQ(ierr);
  ierr = KSPSetInitialGuessNonzero(f->kspU,PETSC_TRUE); CHKERRQ(ierr);
  ierr = KSPSetType(f->kspU,KSPGMRES); CHKERRQ(ierr); //TODO: Try CG?
  ierr = KSPGetPC(f->kspU,&pc); CHKERRQ(ierr);
  ierr = PCSetType(pc,PCSOR); CHKERRQ(ierr);
  PetscReal omega = OptimalSOR(f->da);
  ierr = PCSORSetOmega(pc,omega); CHKERRQ(ierr); //TODO: choose omega by brute force grid technique, verify against algebraic solution
  ierr = PCSORSetIterations(pc,5,1); CHKERRQ(ierr); //TODO: vary number of SOR iterations too, local_its * global_its
  ierr = PCSORSetSymmetric(pc,SOR_LOCAL_SYMMETRIC_SWEEP); CHKERRQ(ierr);
  ierr = KSPSetOptionsPrefix(f->kspU, "vel_"); CHKERRQ(ierr);
  ierr = KSPSetFromOptions(f->kspU); CHKERRQ(ierr); 
  ierr = KSPSetUp(f->kspU); CHKERRQ(ierr);
  
  /* Pressure Solver Setup */
  ierr = KSPCreate(PETSC_COMM_WORLD,&f->kspP); CHKERRQ(ierr);
  ierr = KSPSetOperators(f->kspP,f->matP,f->matP,SAME_PRECONDITIONER); CHKERRQ(ierr);
  ierr = KSPSetInitialGuessNonzero(f->kspP,PETSC_TRUE); CHKERRQ(ierr);
  ierr = KSPDefaultConvergedSetUIRNorm(f->kspP); CHKERRQ(ierr);
  ierr = MatNullSpaceCreate(PETSC_COMM_WORLD,PETSC_TRUE,0,0,&f->nullspace); CHKERRQ(ierr);
  ierr = KSPSetNullSpace(f->kspP,f->nullspace); CHKERRQ(ierr);
  ierr = KSPSetType(f->kspP,KSPGMRES); CHKERRQ(ierr); //TODO: try cg?
  ierr = KSPGetPC(f->kspP,&pc); CHKERRQ(ierr);
  ierr = PCSetType(pc,PCSOR); CHKERRQ(ierr);
  ierr = PCSORSetSymmetric(pc,SOR_LOCAL_SYMMETRIC_SWEEP); CHKERRQ(ierr);
  ierr = PCSORSetOmega(pc,omega); CHKERRQ(ierr); //TODO: choose omega by brute force grid technique, verify against algebraic solution
  ierr = PCSORSetIterations(pc,5,1); CHKERRQ(ierr); //TODO: vary number of SOR iterations too, local_its * global_its
//  ierr = PCFactorSetLevels(pc,5); CHKERRQ(ierr);
//  ierr = PCFactorSetMatOrderingType(pc, MATORDERING_ND); CHKERRQ(ierr);
  ierr = KSPSetOptionsPrefix(f->kspP, "poi_"); CHKERRQ(ierr);
  ierr = KSPSetFromOptions(f->kspP); CHKERRQ(ierr);
  ierr = KSPSetUp(f->kspP); CHKERRQ(ierr);
  
  
  *fluid = f;
  
  PetscLogEventEnd(EVENT_FluidFieldCreate,0,0,0,0);
  PetscFunctionReturn(0);
}

PetscReal OptimalSOR( DA da )
{
  PetscReal rho = 0;
  DALocalInfo info;
  DAGetLocalInfo(da, &info);
  rho = PetscCosScalar( PETSC_PI / info.xm );
  rho+= PetscCosScalar( PETSC_PI / info.ym );
  if( info.zm > 3 ) //3D
  {
    rho+= PetscCosScalar( PETSC_PI / info.zm );
    rho /= 3; //TODO: what is the spectral radius for a 3D problem?
  } else {
    rho /= 2;  
  }

  return 2. / ( 1 + PetscSqrtScalar(1 - rho*rho) ); 
}

#undef __FUNCT__
#define __FUNCT__ "FluidFieldDestroy"
PetscErrorCode FluidFieldDestroy(FluidField f)
{
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  
  ierr = DADestroy(f->da); CHKERRQ(ierr);
  ierr = VecDestroy(f->source); CHKERRQ(ierr);
  ierr = VecDestroy(f->p); CHKERRQ(ierr);
  ierr = VecDestroy(f->px); CHKERRQ(ierr);
  ierr = VecDestroy(f->py); CHKERRQ(ierr);
  ierr = VecDestroy(f->u); CHKERRQ(ierr);
  ierr = VecDestroy(f->v); CHKERRQ(ierr);
  ierr = MatNullSpaceDestroy(f->nullspace); CHKERRQ(ierr);
  ierr = MatDestroy(f->matU); CHKERRQ(ierr);
  ierr = MatDestroy(f->matV); CHKERRQ(ierr);
  ierr = MatDestroy(f->matP); CHKERRQ(ierr);
  ierr = KSPDestroy(f->kspU); CHKERRQ(ierr);
  ierr = KSPDestroy(f->kspP); CHKERRQ(ierr);
  
  if(f->is3D)
  {
    ierr = VecDestroy(f->w); CHKERRQ(ierr);
    ierr = VecDestroy(f->pz); CHKERRQ(ierr);
    ierr = MatDestroy(f->matW); CHKERRQ(ierr);
  }
  
  ierr = PetscFree(f); CHKERRQ(ierr);
  
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "VecWrite"
PetscErrorCode VecWrite( char *str, Vec v )
{
  int s = 256;
  char file[256];
  char temp_dir[256];
  PetscViewer view;
  PetscErrorCode ierr;
  PetscFunctionBegin;
  ierr = PetscGetTmp( PETSC_COMM_WORLD, (char*)&temp_dir, s); CHKERRQ(ierr);
  sprintf(file, "%s/%s.Real64",temp_dir,str); 
  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,file,FILE_MODE_WRITE,&view); CHKERRQ(ierr);
  ierr = VecView(v, view); CHKERRQ(ierr);
  ierr = PetscViewerDestroy(view); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "FluidFieldGetSource"
PetscErrorCode FluidFieldGetSource( FluidField f, Vec *source )
{
  PetscErrorCode ierr;
  
  
  PetscFunctionBegin;
  
  ierr = VecZeroEntries(f->source); CHKERRQ(ierr);
  *source = f->source;
  
  PetscFunctionReturn(0);
}