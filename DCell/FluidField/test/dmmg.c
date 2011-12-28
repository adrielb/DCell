#include "FluidField.h"
#include "petscdmmg.h"

PetscErrorCode ComputeMatrix(DMMG,Mat,Mat);
PetscErrorCode ComputeRHS(DMMG,Vec);

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc, char **args) {

  int nlevels = 3, dof = 1, width = 1, NX = 8;
  DMMG *dmmg;
  DA da;
  PetscErrorCode ierr;
  ierr = DCellInit(); CHKERRQ(ierr);
  ierr = DMMGCreate(PETSC_COMM_WORLD,nlevels,PETSC_NULL,&dmmg);CHKERRQ(ierr);
  ierr = DACreate3d(PETSC_COMM_WORLD,DA_NONPERIODIC,DA_STENCIL_STAR,NX,NX,NX,PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE,dof,width,0,0,0,&da);CHKERRQ(ierr);
  ierr = DMMGSetDM(dmmg,(DM)da);CHKERRQ(ierr);
  ierr = DADestroy(da);CHKERRQ(ierr);

  ierr = DMMGSetKSP(dmmg,ComputeRHS,ComputeMatrix);CHKERRQ(ierr);
//  ierr = DMMGView(dmmg,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
  ierr = DAView( DMMGGetDA(dmmg), PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
  ierr = DMMGSetUp(dmmg);CHKERRQ(ierr);
  ierr = DMMGSolve(dmmg);CHKERRQ(ierr);

  PetscReal norm;
  ierr = MatMult(DMMGGetJ(dmmg),DMMGGetx(dmmg),DMMGGetr(dmmg));CHKERRQ(ierr);
  ierr = VecAXPY(DMMGGetr(dmmg),-1.0,DMMGGetRHS(dmmg));CHKERRQ(ierr);
  ierr = VecNorm(DMMGGetr(dmmg),NORM_2,&norm);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Residual norm %G\n",norm);CHKERRQ(ierr);

  ierr = DMMGDestroy(dmmg);CHKERRQ(ierr);
	ierr = DCellFinalize(); CHKERRQ(ierr);
  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "ComputeRHS"
PetscErrorCode ComputeRHS(DMMG dmmg,Vec b)
{
  PetscErrorCode ierr;
  PetscInt       mx,my,mz;
  PetscScalar    h;

  PetscFunctionBegin;
  ierr = DAGetInfo((DA)dmmg->dm,0,&mx,&my,&mz,0,0,0,0,0,0,0);CHKERRQ(ierr);
  h    = 1.0/((mx-1)*(my-1)*(mz-1));
  ierr = VecSet(b,h);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ComputeMatrix"
PetscErrorCode ComputeMatrix(DMMG dmmg,Mat jac,Mat B)
{
  DA             da = (DA)dmmg->dm;
  PetscErrorCode ierr;
  PetscInt       i,j,k,mx,my,mz,xm,ym,zm,xs,ys,zs;
  PetscScalar    v[7],Hx,Hy,Hz,HxHydHz,HyHzdHx,HxHzdHy;
  MatStencil     row,col[7];

  ierr = DAGetInfo(da,0,&mx,&my,&mz,0,0,0,0,0,0,0);CHKERRQ(ierr);
  Hx = 1.0 / (PetscReal)(mx-1); Hy = 1.0 / (PetscReal)(my-1); Hz = 1.0 / (PetscReal)(mz-1);
  HxHydHz = Hx*Hy/Hz; HxHzdHy = Hx*Hz/Hy; HyHzdHx = Hy*Hz/Hx;
  ierr = DAGetCorners(da,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);

  for (k=zs; k<zs+zm; k++){
    for (j=ys; j<ys+ym; j++){
      for(i=xs; i<xs+xm; i++){
        row.i = i; row.j = j; row.k = k;
  if (i==0 || j==0 || k==0 || i==mx-1 || j==my-1 || k==mz-1){
          v[0] = 2.0*(HxHydHz + HxHzdHy + HyHzdHx);
    ierr = MatSetValuesStencil(B,1,&row,1,&row,v,INSERT_VALUES);CHKERRQ(ierr);
  } else {
    v[0] = -HxHydHz;col[0].i = i; col[0].j = j; col[0].k = k-1;
    v[1] = -HxHzdHy;col[1].i = i; col[1].j = j-1; col[1].k = k;
    v[2] = -HyHzdHx;col[2].i = i-1; col[2].j = j; col[2].k = k;
    v[3] = 2.0*(HxHydHz + HxHzdHy + HyHzdHx);col[3].i = row.i; col[3].j = row.j; col[3].k = row.k;
    v[4] = -HyHzdHx;col[4].i = i+1; col[4].j = j; col[4].k = k;
    v[5] = -HxHzdHy;col[5].i = i; col[5].j = j+1; col[5].k = k;
    v[6] = -HxHydHz;col[6].i = i; col[6].j = j; col[6].k = k+1;
    ierr = MatSetValuesStencil(B,1,&row,7,col,v,INSERT_VALUES);CHKERRQ(ierr);
        }
      }
    }
  }
  ierr = MatAssemblyBegin(B,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(B,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  return 0;
}
