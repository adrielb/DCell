#include "petsc.h"
#include "petscda.h"

PetscErrorCode MatWrite( const char *name, Mat mat );
PetscErrorCode MatMake( int MX, DA *da, Mat *m );


typedef struct {
  PetscReal u,v,w,p;
} Field;

typedef struct {
  PetscReal xx,xy,xz,
               yy,yz,
                  zz;
} StrainTensor;

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc, char **args)
{
  PetscErrorCode  ierr;
  ierr = PetscInitialize(&argc, &args, (char *) 0, ""); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "Start %s\n", __FILE__); CHKERRQ(ierr);
  MPI_Comm comm = PETSC_COMM_WORLD;

  /*
   * DA for velocity [u v w p]
   */

  DA daV;
  Mat mat;
  int MX = 100;

  ierr = PetscPrintf(comm,"Started Assembly\n"); CHKERRQ(ierr);
  ierr = MatMake(MX,&daV, &mat); CHKERRQ(ierr);
  ierr = MatWrite("J",mat); CHKERRQ(ierr);
  ierr = PetscPrintf(comm,"Finished Assembly\n"); CHKERRQ(ierr);

  /*
   * DA for strain tensor
   *  [exx exy exz]
   *  [exy eyy eyz]
   *  [exz eyz ezz]
   */
  DA daE;
  ierr = DACreate3d(comm,DA_NONPERIODIC,DA_STENCIL_STAR,
              -MX,-MX,-MX,
              PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE,
              6,1,0,0,0, &daE); CHKERRQ(ierr);
  Vec e;
  ierr = DACreateGlobalVector(daE,&e); CHKERRQ(ierr);

  Vec rhs, sol;
  ierr = DACreateGlobalVector(daV,&rhs); CHKERRQ(ierr);
  ierr = DACreateGlobalVector(daV,&sol); CHKERRQ(ierr);

  int i,j;
  int xs,ys,xm,ym;
  Field **array;
  ierr = DAVecGetArray(daV,rhs,&array); CHKERRQ(ierr);
  ierr = DAGetCorners(daV,&xs,&ys,PETSC_NULL,&xm,&ym,PETSC_NULL); CHKERRQ(ierr);
  for (j = ys; j < ym; ++j) {
    for (i = xs; i < xm; ++i) {
      if( 3*MX/4 > i && i > MX/4 &&
          3*MX/4 > j && j > MX/4  )
      {
        array[j][i].u = 100.;
        array[j][i].v = 100.;
      }
    }
  }
  ierr = DAVecRestoreArray(daV,rhs,&array); CHKERRQ(ierr);


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

  //Split velocity into component matricies [u], [v], [w]
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
  ierr = PCFieldSplitSetBlockSize(pc,3); CHKERRQ(ierr);
  ierr = PCSetUp(pc); CHKERRQ(ierr);

  /* Set solver for each velocity component
   * Split component velocity as parallel blocks along processors
   * Use direct solver for each block
   */
  int nVel;
  KSP *kspVel, *subksp;
  ierr = PCFieldSplitGetSubKSP(pc,&nVel,&kspVel); CHKERRQ(ierr);
  for( i = 0; i < nVel; i++ ) {
    ierr = KSPSetType(kspVel[i],KSPGMRES); CHKERRQ(ierr);
    ierr = KSPGetPC(kspVel[i],&pc); CHKERRQ(ierr);
    ierr = PCSetType(pc, PCASM); CHKERRQ(ierr);
    ierr = PCASMGetSubKSP(pc,PETSC_NULL,PETSC_NULL,&subksp); CHKERRQ(ierr);
    ierr = KSPSetType(subksp[0],KSPPREONLY); CHKERRQ(ierr);
    ierr = KSPGetPC(subksp[0],&pc); CHKERRQ(ierr);
    ierr = PCSetType(pc,PCCHOLESKY); CHKERRQ(ierr);
    ierr = PCFactorSetMatOrderingType(pc,MATORDERING_ND); CHKERRQ(ierr);
    ierr = PCSetUp(pc); CHKERRQ(ierr);
    ierr = KSPView(kspVel[i],PETSC_VIEWER_DEFAULT); CHKERRQ(ierr);
  }

  printf("\n\n\n\n=============\n\n\n\n");


  for (t = 0; t < tmax; ++t) {
    //AdvectSL each Eij
    //Integrate
  }
  ierr = KSPSolve(ksp,rhs,sol); CHKERRQ(ierr);

  PetscViewer binv;
  ierr = PetscViewerBinaryOpen(comm,"/home/abergman/Research/DCell/temp/uvp.Real64",FILE_MODE_WRITE,&binv); CHKERRQ(ierr);
  ierr = VecView(sol, binv); CHKERRQ(ierr);
  ierr = PetscViewerDestroy(binv); CHKERRQ(ierr);

  ierr = KSPDestroy(ksp); CHKERRQ(ierr);
  ierr = VecDestroy(rhs); CHKERRQ(ierr);
  ierr = VecDestroy(sol); CHKERRQ(ierr);
  ierr = MatDestroy(mat); CHKERRQ(ierr);
  ierr = DADestroy(daV); CHKERRQ(ierr);
  ierr = DADestroy(daE); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "End %s\n", __FILE__); CHKERRQ(ierr);
  ierr = PetscFinalize(); CHKERRQ(ierr);
}

#undef __FUNCT__
#define __FUNCT__ "IntegrateStrainRate"
PetscErrorCode IntegrateStrainRate( DA daV, Vec vecV, DA daE, Vec vecE, PetscReal dh, PetscReal dt )
{
  int i,j,k;
  int xs,ys,zs,
      xm,ym,zm;
  PetscReal dx, dy, dz;
  Field ***vel;
  StrainTensor ***e;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  dx=dy=dz=dh;

  // Update strain rate
  ierr = DAGetCorners(daE,&xs,&ys,&zs,&xm,&ym,&zm); CHKERRQ(ierr);
  ierr = DAVecGetArray(daE,vecE,&e); CHKERRQ(ierr);
  ierr = DAVecGetArray(daV,vecV,&vel); CHKERRQ(ierr);
  for (k = zs; k < zm; ++k) {
    for (j = ys; j < ym; ++j) {
      for (i = xs; i < xm; ++i) {
        e[k][j][i].xx += dt * ( vel[k][j][i].u - vel[k][j][i].u ) / dx;
        e[k][j][i].yy += dt * ( vel[k][j][i].v - vel[k][j][i].v ) / dy;
        e[k][j][i].zz += dt * ( vel[k][j][i].w - vel[k][j][i].w ) / dz;
        e[k][j][i].xy += dt * 0.5 * ( ( vel[k][j][i].u - vel[k][j][i].u ) / dy +
                                      ( vel[k][j][i].v - vel[k][j][i].v ) / dx );
        e[k][j][i].xz += dt * 0.5 * ( ( vel[k][j][i].u - vel[k][j][i].u ) / dz +
                                      ( vel[k][j][i].w - vel[k][j][i].w ) / dx );
        e[k][j][i].yz += dt * 0.5 * ( ( vel[k][j][i].w - vel[k][j][i].w ) / dy +
                                      ( vel[k][j][i].v - vel[k][j][i].v ) / dz );
      }
    }
  }
  ierr = DAVecRestoreArray(daE,vecE,&e); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(daV,vecV,&vel); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ElasticDivergence"
PetscErrorCode ElasticDivergence( DA daE, Vec et, DA daV, Vec rhs, PetscReal dh )
{
  int i,j,k;
  int xs,ys,zs,
      xm,ym,zm;
  PetscReal dx, dy, dz;
  Field ***b;
  StrainTensor ***e;
  Vec etl; //Local strain tensor
  PetscErrorCode ierr;

  PetscFunctionBegin;
  dx=dy=dz=dh;
  ierr = DAGetLocalVector(daE, &etl); CHKERRQ(ierr);
  ierr = DAGlobalToLocalBegin(daE,et,INSERT_VALUES,etl); CHKERRQ(ierr);
  ierr = DAGlobalToLocalEnd(daE,et,INSERT_VALUES,etl); CHKERRQ(ierr);
  ierr = DAVecGetArray(daE,et,&e); CHKERRQ(ierr);
  ierr = DAVecGetArray(daV,rhs,&b); CHKERRQ(ierr);
  ierr = DAGetCorners(daV,&xs,&ys,&zs,&xm,&ym,&zm); CHKERRQ(ierr);
  for (k = zs; k < zm; ++k) {
    for (j = ys; j < ym; ++j) {
      for (i = xs; i < xm; ++i) {
        // Del.e
        b[k][j][i].u =  ( e[k][j][i].xx - e[k][j][i].xx ) / dx +
                        ( e[k][j][i].xy - e[k][j][i].xy ) / dy +
                        ( e[k][j][i].xz - e[k][j][i].xz ) / dz;
        b[k][j][i].v =  ( e[k][j][i].xy - e[k][j][i].xy ) / dx +
                        ( e[k][j][i].yy - e[k][j][i].yy ) / dy +
                        ( e[k][j][i].yz - e[k][j][i].yz ) / dz;
        b[k][j][i].w =  ( e[k][j][i].xz - e[k][j][i].xz ) / dx +
                        ( e[k][j][i].yz - e[k][j][i].yz ) / dy +
                        ( e[k][j][i].zz - e[k][j][i].zz ) / dz;
      }
    }
  }
  ierr = DAVecRestoreArray(daV,rhs,&b); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(daE,et,&e); CHKERRQ(ierr);
  ierr = DARestoreLocalVector(daE, &etl); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MatMake"
PetscErrorCode MatMake( int MX, DA *d, Mat *m )
{
  Mat mat;
  DA da;
  const int dof = 4; // [u v w p]
  const int P = 3; // index of pressure component
  DALocalInfo i;
  DAGetLocalInfo(da, &i);
  PetscReal dd = 1./MX;
  PetscReal hx = dd, hy = dd, hz = dd;
  PetscReal dh[3] = {hx,hy,hz};
  const int n=7; // Length of col[7];
  int c; //dummy index for looping over cols[]
  MatStencil row, col[n], pcol[2];
  PetscReal pval[2],
            val[] = {-hx*hx,-hx*hx,
                      -hy*hy,-hy*hy,
                      -hz*hz,-hz*hz, 2*(hx*hx + hy*hy + hz*hz)},
            val_D[] = { 0, 0, 0,
                        0, 0, 0, 1},
            val_div[6] = {hx,-hx,
                          hy,-hy,
                          hz,-hz};
  int end[][3] = {{1,2,2}, // u
                  {2,1,2}, // v
                  {2,2,1}};// w
  int op[][3] = {{1,0,0},
                 {0,1,0},
                 {0,0,1}};
  PetscErrorCode ierr;

  PetscFunctionBegin;

  ierr = DACreate3d(PETSC_COMM_WORLD,DA_NONPERIODIC,DA_STENCIL_STAR,
            -MX,-MX,-MX,
            PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE,
            dof,1,
            0,0,0, &da); CHKERRQ(ierr);

  ierr = DAGetMatrix(da, MATMPIAIJ, &mat); CHKERRQ(ierr);

  pcol[0].c = P;
  pcol[1].c = P;

  for(     row.k = i.zs; row.k < i.zs+i.zm; ++row.k) {
    for(   row.j = i.ys; row.j < i.ys+i.ym; ++row.j) {
      for( row.i = i.xs; row.i < i.xs+i.xm; ++row.i) {
        row.c = 0;
        for( c = 0; c < n; ++c)
          PetscMemcpy(&col[c],&row,sizeof(MatStencil));
        col[0].i--; col[1].i++;  // [-hx, hx]
        col[2].j--; col[3].j++;  // [-hy, hy]
        col[4].k--; col[5].k++;  // [-hz, hz]
        if( col[1].i >= i.mx ) col[1].i = -1;
        if( col[3].j >= i.my ) col[3].j = -1;
        if( col[5].k >= i.mz ) col[5].k = -1;

        // Laplace([u,v,w]) - grad([px,py,pz])
        for( row.c = 0; row.c < 3; row.c++ )  // [u,v,w]
        {
          for( c = 0; c < n; ++c)
            col[c].c = row.c;

          if( 0 < row.i && row.i < i.mx - end[row.c][0] &&
              0 < row.j && row.j < i.my - end[row.c][1] &&
              0 < row.k && row.k < i.mz - end[row.c][2] ) {

            // Laplace(u)
            ierr = MatSetValuesStencil(mat,1,&row,n,col,val,INSERT_VALUES); CHKERRQ(ierr);

            // px[i,j,k] = ( p[i,j,k] - p[i-1,j,k] ) / hx
            pval[0] = -dh[row.c];
            pval[1] =  dh[row.c];

            pcol[0].i = row.i;
            pcol[0].j = row.j;
            pcol[0].k = row.k;
            pcol[1].i = row.i - op[row.c][0];
            pcol[1].j = row.j - op[row.c][1];
            pcol[1].k = row.k - op[row.c][2];
            ierr = MatSetValuesStencil(mat,1,&row,2,pcol,pval,INSERT_VALUES); CHKERRQ(ierr);
          } else {
            ierr = MatSetValuesStencil(mat,1,&row,n,col,val_D,INSERT_VALUES); CHKERRQ(ierr);
          }
        }

        row.c = P;
        for( c = 0; c < n; ++c)
          col[c].c = row.c;

        if( row.j == 0 && (row.i == 0 || row.i == i.mx-2) ) {
          //Neumann BC on lower corners
          pval[0] =  hy;
          pval[1] = -hy;

          pcol[0].i = row.i;
          pcol[0].j = row.j+1;
          pcol[0].k = row.k;
          pcol[1].i = row.i;
          pcol[1].j = row.j;
          pcol[1].k = row.k;
          ierr = MatSetValuesStencil(mat,1,&row,2,pcol,pval,INSERT_VALUES); CHKERRQ(ierr);
        } else if( row.j == i.mx - 2 && (row.i == 0 || row.i == i.my-2) ) {
          //Neumann BC on upper corners
          pval[0] =  hy;
          pval[1] = -hy;

          pcol[0].i = row.i;
          pcol[0].j = row.j;
          pcol[0].k = row.k;
          pcol[1].i = row.i;
          pcol[1].j = row.j-1;
          pcol[1].k = row.k;
          ierr = MatSetValuesStencil(mat,1,&row,2,pcol,pval,INSERT_VALUES); CHKERRQ(ierr);
        } else if( 0 <= row.i && row.i < i.mx-1 &&
                   0 <= row.j && row.j < i.my-1 &&
                   0 <= row.k && row.k < i.mz-1 ) {
          //Divergence equation for interior nodes
          col[0].c = 0;
          col[0].i = row.i+1;
          col[1].c = 0;
          col[1].i = row.i;
          col[2].c = 1;
          col[2].j = row.j+1;
          col[3].c = 1;
          col[3].j = row.j;
          col[4].c = 2;
          col[4].k = row.k+1;
          col[5].c = 2;
          col[5].k = row.k;
          ierr = MatSetValuesStencil(mat,1,&row,6,col,val_div,INSERT_VALUES); CHKERRQ(ierr);
        } else if( row.i == i.mx-1 || row.j == i.my-1 ) {
          //Overhanging pressure nodes on upper boundaries
          ierr = MatSetValuesStencil(mat,1,&row,n,col,val_D,INSERT_VALUES); CHKERRQ(ierr);
        }
      }
    }
  }

  ierr = MatAssemblyBegin(mat,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(mat,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
  /*
   *  Make velocity block symmetric by eliminating BCs
   */
  for(     row.k = i.zs; row.k < i.zs+i.zm; ++row.k) {
    for(   row.j = i.ys; row.j < i.ys+i.ym; ++row.j) {
      for( row.i = i.xs; row.i < i.xs+i.xm; ++row.i) {
        row.c = 0;
        for( c = 0; c < n; ++c)
          PetscMemcpy(&col[c],&row,sizeof(MatStencil));
        col[0].i--; col[1].i++;  // [-hx, hx]
        col[2].j--; col[3].j++;  // [-hy, hy]
        col[4].k--; col[5].k++;  // [-hz, hz]
        if( col[1].i >= i.mx ) col[1].i = -1;
        if( col[3].j >= i.my ) col[3].j = -1;
        if( col[5].k >= i.mz ) col[5].k = -1;

        for( row.c = 0; row.c < 3; row.c++ )  // [u,v,w]
        {
          for( c = 0; c < n; ++c)
            col[c].c = row.c;

          if( 0 < row.i && row.i < i.mx - end[row.c][0] &&
              0 < row.j && row.j < i.my - end[row.c][1] &&
              0 < row.k && row.k < i.mz - end[row.c][2] ) {
          } else {
            ierr = MatSetValuesStencil(mat,n,col,1,&row,val_D,INSERT_VALUES); CHKERRQ(ierr);
          }
        } // row.c
      } // row.i
    } // row.j
  } // row.k
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
