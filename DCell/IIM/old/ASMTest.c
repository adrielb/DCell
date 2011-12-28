#include "petscksp.h"
#include "petscda.h"
#include "Main.h"

PetscErrorCode viz_apply( void *ctx, Vec in, Vec out);

PetscErrorCode PetscMain()
{
  DA da;
  Mat A;
  DALocalInfo i;
  MatStencil row, col[5];
  PetscReal val[5] = {1,1,1,1,-5};
  int c, rank;
  Vec sol;
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
  val[4]=rank;
  ierr = DACreate2d(PETSC_COMM_WORLD,DA_NONPERIODIC,DA_STENCIL_BOX,
      -1,-1,PETSC_DECIDE,PETSC_DECIDE,1,1,0,0, &da); CHKERRQ(ierr);
  ierr = DAGetMatrix(da, MATMPIAIJ,&A); CHKERRQ(ierr);
  ierr = DAGetLocalInfo(da, &i); CHKERRQ(ierr);
  
  
  for( row.j = i.ys; row.j < i.ys+i.ym; ++row.j)
  {
    for( row.i = i.xs; row.i < i.xs+i.xm; ++row.i)
    {
      for( c = 0; c < 5; ++c)
      {
        PetscMemcpy(&col[c],&row,sizeof(MatStencil));
      }
      col[0].i--; col[2].j--; 
      col[1].i++; col[3].j++;
      if( col[1].i == i.mx ) col[1].i = -1;
      if( col[3].j == i.my ) col[3].j = -1;
//      val[4] = row.i + 0.1*row.j;
      val[4] = rank;
      ierr = MatSetValuesStencil(A,1,&row,5,col,val,INSERT_VALUES); CHKERRQ(ierr);
    }
  }
  
  ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  
  PetscViewer view;
  PetscViewerASCIIOpen(PETSC_COMM_SELF, "../temp/mat.dat", &view);
  MatView(A,view);
  PetscViewerDestroy(view);
  
  ierr = DACreateGlobalVector(da, &sol); CHKERRQ(ierr);
  Vec local;
  ierr = DAGetLocalVector(da, &local);; CHKERRQ(ierr);
  
  int k,j;
  PetscReal **p;
  ierr = DAVecGetArray(da,local,&p); CHKERRQ(ierr);
  for( j = i.ys; j < i.ys+i.ym; ++j)
  {
    for( k = i.xs; k < i.xs+i.xm; ++k)
    {
      p[j][k] = k + 0.1 * j;
    }
  }
  ierr = VecSet(local,rank); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(da,local,&p); CHKERRQ(ierr);
  ierr = DALocalToGlobal(da,local,INSERT_VALUES,sol); CHKERRQ(ierr);
  ierr = DARestoreLocalVector(da,&local); CHKERRQ(ierr);
  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD, "../temp/vec.dat", FILE_MODE_WRITE, &view);; CHKERRQ(ierr);
  ierr = VecView(sol,view);; CHKERRQ(ierr);
  ierr = PetscViewerDestroy(view); CHKERRQ(ierr);

  KSP ksp, *subksp;
  PC pc, subpc;
  ierr = KSPCreate(PETSC_COMM_WORLD,&ksp); CHKERRQ(ierr);
  ierr = KSPSetOperators(ksp,A,A,DIFFERENT_NONZERO_PATTERN); CHKERRQ(ierr);
  ierr = KSPGetPC(ksp,&pc); CHKERRQ(ierr);
  ierr = PCSetType(pc,PCASM); CHKERRQ(ierr);
  ierr = PCASMSetOverlap(pc,1); CHKERRQ(ierr);
  ierr = KSPSetUp(ksp); CHKERRQ(ierr);
  ierr = PCASMGetSubKSP(pc,0,0,&subksp); CHKERRQ(ierr);
  ierr = KSPSetType(subksp[0],KSPPREONLY); CHKERRQ(ierr);
  ierr = KSPGetPC(subksp[0],&subpc); CHKERRQ(ierr);
  ierr = PCSetType(subpc,PCSHELL); CHKERRQ(ierr);
  ierr = PCShellSetApply(subpc,viz_apply); CHKERRQ(ierr);

  int overlap, xstart, xend, ystart, yend;
  overlap = 1;
  xstart = i.xs - overlap;
  xstart = xstart < 0 ? 0 : xstart;
  xend = i.xs+i.xm + overlap;
  xend = xend <= i.mx ? xend : i.mx;
  ystart = i.ys - overlap;
  ystart = ystart < 0 ? 0 : ystart;
  yend = i.ys+i.ym + overlap;
  yend = yend <= i.my ? yend : i.my;
  int idx_size = (i.xm + 2*overlap) * (i.ym + 2*overlap);
  int *idx;
  int count = 0;
  PetscMalloc(idx_size*sizeof(int),&idx);

  for( j = ystart; j < yend; ++j)
  {
    for( k = xstart; k < xend; ++k)
    {
      idx[count] = k + i.mx * j;
      count++;
    }
  }
  PetscPrintf(PETSC_COMM_SELF, "[%d]xs: %d\txe: %d\tys: %d\tye: %d\tidx_size: %d\tcount: %d\n",rank,xstart,xend,ystart,yend,idx_size,count);
  
  IS subIS;
  ierr = ISCreateGeneralNC(PETSC_COMM_SELF,count,idx,&subIS); CHKERRQ(ierr);
//  ierr = PCASMSetLocalSubdomains(pc,1,&subIS); CHKERRQ(ierr);
  
  
  PetscInt n_local;
  Mat *mats, pmat;
  ierr = PCASMGetLocalSubmatrices(pc,&n_local,&mats); CHKERRQ(ierr);
  ierr = MatPermute(mats[0],subIS,subIS,&pmat); CHKERRQ(ierr);
//  printf("N_LOCAL: %d\n",n_local);
  char mat_file[128];
  sprintf(mat_file,"../temp/mat.%d.dat",rank);
  PetscViewerASCIIOpen(PETSC_COMM_SELF, mat_file, &view);
  MatView(pmat,view);
  PetscViewerDestroy(view);
  
  /*
   * Destroy's distributed matrix
   * 
  KSPSetType(subksp,KSPPREONLY);
  PCASMSetUseInPlace(pc)
  */
  
  ierr = KSPDestroy(ksp); CHKERRQ(ierr);
  ierr = MatDestroy(A); CHKERRQ(ierr);
  ierr = DADestroy(da); CHKERRQ(ierr);
  ierr = VecDestroy(sol); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


PetscErrorCode viz_apply( void *ctx, Vec in, Vec out)
{
  PetscViewer view;
  PetscErrorCode ierr;
  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD, "../temp/vec.dat", FILE_MODE_WRITE, &view); CHKERRQ(ierr);
  ierr = VecView(in,view); CHKERRQ(ierr);
  ierr = PetscViewerDestroy(view); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

int RunCheck()
{
  return 0;
}