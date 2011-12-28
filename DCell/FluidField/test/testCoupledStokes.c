#include "petsc.h"
#include "petscda.h"

PetscErrorCode MatWrite( const char *name, Mat mat );
PetscErrorCode MatMake( int MX, int dof, DA *da, Mat *m );

typedef struct {
  PetscReal u,v,p;
} Field;

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc, char **args)
{
  PetscErrorCode  ierr;
  ierr = PetscInitialize(&argc, &args, (char *) 0, ""); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "Start %s\n", __FILE__); CHKERRQ(ierr);
  MPI_Comm comm = PETSC_COMM_WORLD;

  DA da;
  Mat mat;
  int MX = 100;
  int dof = 3;

  ierr = PetscPrintf(comm,"Started Assembly\n"); CHKERRQ(ierr);
  ierr = MatMake(MX,dof,&da, &mat); CHKERRQ(ierr);
  ierr = MatWrite("J",mat); CHKERRQ(ierr);
  ierr = PetscPrintf(comm,"Finished Assembly\n"); CHKERRQ(ierr);

  Vec rhs, sol;
  ierr = DACreateGlobalVector(da,&rhs); CHKERRQ(ierr);
  ierr = DACreateGlobalVector(da,&sol); CHKERRQ(ierr);

  int i,j;
  int xs,ys,xm,ym;
  Field **array;
  ierr = DAVecGetArray(da,rhs,&array); CHKERRQ(ierr);
  ierr = DAGetCorners(da,&xs,&ys,PETSC_NULL,&xm,&ym,PETSC_NULL); CHKERRQ(ierr);
  for (j = xs; j < xm; ++j) {
    for (i = ys; i < ym; ++i) {
      if( 3*MX/4 > i && i > MX/4 &&
          3*MX/4 > j && j > MX/4  )
      {
        array[j][i].u = 100.;
        array[j][i].v = 100.;
      }
    }
  }
  ierr = DAVecRestoreArray(da,rhs,&array); CHKERRQ(ierr);


  KSP ksp;
  ierr = KSPCreate(comm,&ksp); CHKERRQ(ierr);
  ierr = KSPSetType(ksp,KSPPREONLY); CHKERRQ(ierr);
  ierr = KSPSetOperators(ksp,mat,mat,DIFFERENT_NONZERO_PATTERN); CHKERRQ(ierr);
  ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);

  //Split pressure from velocity
  PC pc;
  ierr = KSPGetPC(ksp,&pc); CHKERRQ(ierr);
  ierr = PCSetType(pc, PCFIELDSPLIT); CHKERRQ(ierr);
  ierr = PCFieldSplitSetType(pc,PC_COMPOSITE_SCHUR); CHKERRQ(ierr);
  ierr = PCFieldSplitSetBlockSize(pc,3); CHKERRQ(ierr);
  ierr = PCFieldSplitSetFields(pc,2,(int[]){0,1}); CHKERRQ(ierr);
  ierr = PCFieldSplitSetFields(pc,1,(int[]){2}); CHKERRQ(ierr);
  ierr = PCSetUp(pc); CHKERRQ(ierr);

  //Split velocity into u-v component matricies
  int nVelP;
  KSP *kspVelP;
  ierr = PCFieldSplitGetSubKSP(pc,&nVelP,&kspVelP); CHKERRQ(ierr);
  ierr = KSPSetTolerances(kspVelP[1],1e-50,1e-3,PETSC_DEFAULT,PETSC_DEFAULT); CHKERRQ(ierr);
//  ierr = KSPSetType(kspVelP[1],KSPCG); CHKERRQ(ierr);
  ierr = KSPGetPC(kspVelP[1],&pc); CHKERRQ(ierr);
  ierr = PCSetType(pc,PCNONE); CHKERRQ(ierr);
  ierr = KSPSetType(kspVelP[0],KSPPREONLY); CHKERRQ(ierr);
  ierr = KSPGetPC(kspVelP[0],&pc); CHKERRQ(ierr);
  ierr = PCSetType(pc, PCFIELDSPLIT); CHKERRQ(ierr);
  ierr = PCFieldSplitSetType(pc,PC_COMPOSITE_ADDITIVE); CHKERRQ(ierr);
  ierr = PCFieldSplitSetBlockSize(pc,2); CHKERRQ(ierr);
  ierr = PCSetUp(pc); CHKERRQ(ierr);

  //Set solver for each velocity component
  int nVel;
  KSP *kspVel;
  ierr = PCFieldSplitGetSubKSP(pc,&nVel,&kspVel); CHKERRQ(ierr);
  for( i = 0; i < nVel; i++ ) {
    // Direct Method
    ierr = KSPSetType(kspVel[i],KSPPREONLY); CHKERRQ(ierr);
    ierr = KSPGetPC(kspVel[i],&pc); CHKERRQ(ierr);
    ierr = PCSetType(pc,PCCHOLESKY); CHKERRQ(ierr);
    ierr = PCFactorSetMatOrderingType(pc,MATORDERING_ND); CHKERRQ(ierr);
    ierr = PCSetUp(pc); CHKERRQ(ierr);
    //
    /* Iterative Method
    ierr = KSPSetType(kspVel[i],KSPCG); CHKERRQ(ierr);
    ierr = KSPGetPC(kspVel[i],&pc); CHKERRQ(ierr);
//    ierr = KSPSetTolerances(kspVel[i],1e-3,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT); CHKERRQ(ierr);
    ierr = PCSetType(pc, PCICC); CHKERRQ(ierr);
    ierr = PCFactorSetLevels(pc, 3); CHKERRQ(ierr);
    ierr = PCFactorSetMatOrderingType(pc,MATORDERING_ND); CHKERRQ(ierr);
    ierr = KSPView(kspVel[i],PETSC_VIEWER_DEFAULT); CHKERRQ(ierr);
    */
  }

//  ierr = KSPView(ksp,PETSC_VIEWER_DEFAULT); CHKERRQ(ierr);

  //Split component velocity as parallel blocks along processors


  /*
  Mat M,G,D,Z;
  ierr = PCFieldSplitGetSchurBlocks(pc,&M,&G,&D,&Z); CHKERRQ(ierr);
  ierr = MatWrite("M",M); CHKERRQ(ierr);
  ierr = MatWrite("G",G); CHKERRQ(ierr);
  ierr = MatWrite("D",D); CHKERRQ(ierr);
  ierr = MatWrite("Z",Z); CHKERRQ(ierr);
  */

  printf("\n\n\n\n=============\n\n\n\n");

  ierr = KSPSolve(ksp,rhs,sol); CHKERRQ(ierr);

  PetscViewer binv;
  ierr = PetscViewerBinaryOpen(comm,"/home/abergman/Research/DCell/temp/uvp.Real64",FILE_MODE_WRITE,&binv); CHKERRQ(ierr);
  ierr = VecView(sol, binv); CHKERRQ(ierr);
  ierr = PetscViewerDestroy(binv); CHKERRQ(ierr);

  ierr = KSPDestroy(ksp); CHKERRQ(ierr);
  ierr = VecDestroy(rhs); CHKERRQ(ierr);
  ierr = VecDestroy(sol); CHKERRQ(ierr);
  ierr = MatDestroy(mat); CHKERRQ(ierr);
  ierr = DADestroy(da); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "End %s\n", __FILE__); CHKERRQ(ierr);
  ierr = PetscFinalize(); CHKERRQ(ierr);
}

#undef __FUNCT__
#define __FUNCT__ "MatMake"
PetscErrorCode MatMake( int MX, int dof, DA *d, Mat *m )
{
  Mat mat;
  DA da;
  PetscErrorCode ierr;

  PetscFunctionBegin;

  ierr = DACreate2d(PETSC_COMM_WORLD,DA_NONPERIODIC,DA_STENCIL_STAR,
            -MX,-MX,PETSC_DECIDE,PETSC_DECIDE,dof,1,0,0, &da); CHKERRQ(ierr);

//  ierr = DAGetMatrix(da, MATMPIAIJ, &mat); CHKERRQ(ierr);
  ierr = DAGetMatrix(da, MATSEQAIJ, &mat); CHKERRQ(ierr);

  DALocalInfo i;
  DAGetLocalInfo(da, &i);
  PetscReal hx = 1./MX, hy = 1./MX;
  MatStencil row, col[5], pcol[2];
  PetscReal pval[2],
            val[5] = {-hx*hx,-hx*hx,-hy*hy,-hy*hy,2*(hx*hx+hy*hy)},
            val_D[5] = { 0, 0, 0, 0, 1},
            val_div[4] = {hx,-hx,hy,-hy};
  int c;
  for( row.j = i.ys; row.j < i.ys+i.ym; ++row.j)
  {
    for( row.i = i.xs; row.i < i.xs+i.xm; ++row.i)
    {
      row.c = 0;
      for( c = 0; c < 5; ++c)
        PetscMemcpy(&col[c],&row,sizeof(MatStencil));
      col[0].i--; col[1].i++;  // [-hx, hx]
      col[2].j--; col[3].j++;  // [-hy, hy]
      if( col[1].i >= i.mx ) col[1].i = -1;
      if( col[3].j >= i.my ) col[3].j = -1;

      // Laplace(u) - px
      if( 0 < row.i && row.i < i.mx - 1 &&
          0 < row.j && row.j < i.my - 2 ) {

        ierr = MatSetValuesStencil(mat,1,&row,5,col,val,INSERT_VALUES); CHKERRQ(ierr);

        // px[i,j] = ( p[i+1,j] - p[i,j] ) / hx
        pval[0] = -hx;
        pval[1] =  hx;

        pcol[0].c = 2;
        pcol[0].i = row.i;
        pcol[0].j = row.j;
        pcol[1].c = 2;
        pcol[1].i = row.i-1;
        pcol[1].j = row.j;
        ierr = MatSetValuesStencil(mat,1,&row,2,pcol,pval,INSERT_VALUES); CHKERRQ(ierr);
      } else {
        ierr = MatSetValuesStencil(mat,1,&row,5,col,val_D,INSERT_VALUES); CHKERRQ(ierr);
      }

      row.c = 1;
      for( c = 0; c < 5; ++c)
        col[c].c = row.c;
      // Laplace(v) - py
      if( 0 < row.j && row.j < i.my - 1 &&
          0 < row.i && row.i < i.mx - 2 ) {

        ierr = MatSetValuesStencil(mat,1,&row,5,col,val,INSERT_VALUES); CHKERRQ(ierr);

        // py[i,j] = ( p[i,j+1] - p[i,j] ) / hy
        pval[0] = -hy;
        pval[1] =  hy;

        pcol[0].c = 2;
        pcol[0].i = row.i;
        pcol[0].j = row.j;
        pcol[1].c = 2;
        pcol[1].i = row.i;
        pcol[1].j = row.j-1;
        ierr = MatSetValuesStencil(mat,1,&row,2,pcol,pval,INSERT_VALUES); CHKERRQ(ierr);
      } else {
        ierr = MatSetValuesStencil(mat,1,&row,5,col,val_D,INSERT_VALUES); CHKERRQ(ierr);
      }

      row.c = 2;
      for( c = 0; c < 5; ++c)
        col[c].c = row.c;

      if( row.j == 0 && (row.i == 0 || row.i == i.mx-2) ) {
        //Neumann BC on lower corners
        pval[0] =  hy;
        pval[1] = -hy;

        pcol[0].c = 2;
        pcol[0].i = row.i;
        pcol[0].j = row.j+1;
        pcol[1].c = 2;
        pcol[1].i = row.i;
        pcol[1].j = row.j;
        ierr = MatSetValuesStencil(mat,1,&row,2,pcol,pval,INSERT_VALUES); CHKERRQ(ierr);
      } else if( row.j == i.mx - 2 && (row.i == 0 || row.i == i.my-2) ) {
        //Neumann BC on upper corners
        pval[0] =  hy;
        pval[1] = -hy;

        pcol[0].c = 2;
        pcol[0].i = row.i;
        pcol[0].j = row.j;
        pcol[1].c = 2;
        pcol[1].i = row.i;
        pcol[1].j = row.j-1;
        ierr = MatSetValuesStencil(mat,1,&row,2,pcol,pval,INSERT_VALUES); CHKERRQ(ierr);
      } else if( 0 <= row.i && row.i < i.mx-1 &&
                 0 <= row.j && row.j < i.my-1 ) {
        //Divergence equation for interior nodes
        col[0].c = 0;
        col[0].i = row.i+1;
        col[1].c = 0;
        col[1].i = row.i;
        col[2].c = 1;
        col[2].j = row.j+1;
        col[3].c = 1;
        col[3].j = row.j;
        ierr = MatSetValuesStencil(mat,1,&row,4,col,val_div,INSERT_VALUES); CHKERRQ(ierr);
      } else if( row.i == i.mx-1 || row.j == i.my-1 ) {
        //Overhanging pressure nodes on upper boundaries
        ierr = MatSetValuesStencil(mat,1,&row,5,col,val_D,INSERT_VALUES); CHKERRQ(ierr);
      }
    }
  }

  ierr = MatAssemblyBegin(mat,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(mat,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
  /*
   *  Make velocity block symmetric by eliminating BCs
   */
  for( row.j = i.ys; row.j < i.ys+i.ym; ++row.j)
  {
    for( row.i = i.xs; row.i < i.xs+i.xm; ++row.i)
    {
      row.c = 0;
      for( c = 0; c < 5; ++c)
        PetscMemcpy(&col[c],&row,sizeof(MatStencil));
      col[0].i--; col[1].i++;  // [-hx, hx]
      col[2].j--; col[3].j++;  // [-hy, hy]
      if( col[1].i >= i.mx ) col[1].i = -1;
      if( col[3].j >= i.my ) col[3].j = -1;

      // Laplace(u) - px
      if( 0 < row.i && row.i < i.mx - 1 &&
          0 < row.j && row.j < i.my - 2 ) {

      } else {
        ierr = MatSetValuesStencil(mat,5,col,1,&row,val_D,INSERT_VALUES); CHKERRQ(ierr);
      }

      row.c = 1;
      for( c = 0; c < 5; ++c)
        col[c].c = row.c;
      // Laplace(v) - py
      if( 0 < row.j && row.j < i.my - 1 &&
          0 < row.i && row.i < i.mx - 2 ) {

      } else {
        ierr = MatSetValuesStencil(mat,5,col,1,&row,val_D,INSERT_VALUES); CHKERRQ(ierr);
      }
    }
  }

  ierr = MatAssemblyBegin(mat,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(mat,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

  *m = mat;
  *d = da;

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MatWrite"
PetscErrorCode MatWrite( const char *name, Mat mat )
{
  PetscViewer view;
  int dir_len = 512;
  char tempdir[512], dir[512];
  PetscErrorCode ierr;

  PetscFunctionBegin;

  ierr = PetscGetTmp( PETSC_COMM_WORLD, tempdir, dir_len); CHKERRQ(ierr);
  ierr = PetscSNPrintf(dir, dir_len, "%s/%s.mat", tempdir, name); CHKERRQ(ierr);
  ierr = PetscViewerASCIIOpen(PETSC_COMM_SELF, dir, &view); CHKERRQ(ierr);
  ierr = MatView(mat, view); CHKERRQ(ierr);
  ierr = PetscViewerDestroy(view); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}
