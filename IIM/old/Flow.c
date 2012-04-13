#include "Flow.h"
#include "mkl.h"
#include "LevelSet.h"

#undef __FUNCT__
#define __FUNCT__ "PetscMain"
PetscErrorCode PetscMain(  )
{
	PetscErrorCode ierr;
	DMMG *dmmg;
  DA da;
  PC pc, pcmg;
  KSP ksp;
	
	PetscFunctionBegin;
	PetscLogEventRegister(&EVENT_ComputeRHS,"ComputeRHS", 0);
	PetscLogEventRegister(&EVENT_ComputeJacobian,"ComputeJacobian", 0);
	
  DMMGCreate(PETSC_COMM_SELF, 3, PETSC_NULL, &dmmg);
  ierr = DACreate3d(PETSC_COMM_SELF,//MPI Communicator   
    DA_NONPERIODIC,         //DA_NONPERIODIC, DA_XPERIODIC, DA_YPERIODIC, DA_XYPERIODIC
    DA_STENCIL_STAR,       //DA_STENCIL_BOX or DA_STENCIL_STAR
    32, 32, 32,    //Global dimension
    1, 1, 1,   //Number procs per dim
    1,        //dof
    1,       //stencil width
    0, 0, 0,//specific array of nodes
    &da); CHKERRQ(ierr);
  
  ierr = DMMGSetDM( dmmg, (DM) da ); CHKERRQ(ierr);
  ierr = DADestroy(da); CHKERRQ(ierr);
  
  ierr = DMMGSetKSP(dmmg, ComputeRHS, ComputeJacobian); CHKERRQ(ierr);
  ierr = DMMGSetUp(dmmg); CHKERRQ(ierr);
  
  ksp = DMMGGetKSP(dmmg);
  KSPSetType(ksp, KSPPREONLY);
  KSPGetPC(ksp, &pcmg);
  for( int i = 0; i < DMMGGetLevels(dmmg); i++)
  {
  	ierr = PCMGGetSmoother(pcmg, i, &ksp); CHKERRQ(ierr);
  	KSPSetType(ksp, KSPCG);
		KSPGetPC(ksp, &pc);
		PCSetType(pc, PCICC);
//		PCFactorSetLevels(pc, 3);
  }
  ierr = DMMGSolve(dmmg); CHKERRQ(ierr);
  ierr = DMMGView(dmmg, PETSC_VIEWER_STDOUT_SELF); CHKERRQ(ierr);
  
  DMMGDestroy(dmmg); 
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ComputeRHS"
PetscErrorCode ComputeRHS( DMMG dmmg,Vec b )
{
	PetscErrorCode ierr;
	PetscInt       mx,my,mz;
	PetscScalar    h;
	
	PetscFunctionBegin;
	PetscLogEventBegin(EVENT_ComputeRHS,0,0,0,0);
	
  DAGetInfo((DA)dmmg->dm,0,&mx,&my,&mz,0,0,0,0,0,0,0);
  h = 1.0/((mx-1)*(my-1)*(mz-1));
  VecSet(b,h);
		
	PetscLogEventEnd(EVENT_ComputeRHS,0,0,0,0);
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ComputeJacobian"
PetscErrorCode ComputeJacobian( DMMG dmmg,Mat jac,Mat B )
{
	PetscErrorCode ierr;
	DA             da = (DA)dmmg->dm;
  PetscInt       i,j,k,mx,my,mz,xm,ym,zm,xs,ys,zs;
  PetscScalar    v[7],Hx,Hy,Hz,HxHydHz,HyHzdHx,HxHzdHy;
  MatStencil     row,col[7];
	
	PetscFunctionBegin;
	PetscLogEventBegin(EVENT_ComputeJacobian,0,0,0,0);
	
	DAGetInfo(da,0,&mx,&my,&mz,0,0,0,0,0,0,0);
  Hx = 1.0 / (PetscReal)(mx-1); Hy = 1.0 / (PetscReal)(my-1); Hz = 1.0 / (PetscReal)(mz-1);
  HxHydHz = Hx*Hy/Hz; HxHzdHy = Hx*Hz/Hy; HyHzdHx = Hy*Hz/Hx;
  DAGetCorners(da,&xs,&ys,&zs,&xm,&ym,&zm);
   
  for (k=zs; k<zs+zm; k++){
    for (j=ys; j<ys+ym; j++){
      for(i=xs; i<xs+xm; i++){
        row.i = i; row.j = j; row.k = k;
        if (i==0 || j==0 || k==0 || i==mx-1 || j==my-1 || k==mz-1)
        {
          v[0] = 2.0*(HxHydHz + HxHzdHy + HyHzdHx);
          MatSetValuesStencil(B,1,&row,1,&row,v,INSERT_VALUES);
        } else {
          v[0] = -HxHydHz;col[0].i = i; col[0].j = j; col[0].k = k-1;
	        v[1] = -HxHzdHy;col[1].i = i; col[1].j = j-1; col[1].k = k;
  	      v[2] = -HyHzdHx;col[2].i = i-1; col[2].j = j; col[2].k = k;
    	    v[3] = 2.0*(HxHydHz + HxHzdHy + HyHzdHx);col[3].i = row.i; col[3].j = row.j; col[3].k = row.k;
      	  v[4] = -HyHzdHx;col[4].i = i+1; col[4].j = j; col[4].k = k;
        	v[5] = -HxHzdHy;col[5].i = i; col[5].j = j+1; col[5].k = k;
        	v[6] = -HxHydHz;col[6].i = i; col[6].j = j; col[6].k = k+1;
        	MatSetValuesStencil(B,1,&row,7,col,v,INSERT_VALUES);
        }
      }
    }
  }
  MatAssemblyBegin(B,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(B,MAT_FINAL_ASSEMBLY);
//	MatView(B, PETSC_VIEWER_STDOUT_SELF);
	PetscLogEventEnd(EVENT_ComputeJacobian,0,0,0,0);
	PetscFunctionReturn(0);
}



double SurfaceDerivative()

// Equ 1.34
LocalCoor( PetscReal X, PetscReal Y)
{
  c = cos(theta);
  s = sin(theta);
  dx = x - X;
  dy = y - Y;
  
  xi =    dx * c + dy * s
  eta=-1.*dx * s + dy * c  
}

int RunCheck()
{
	return 0;
}
