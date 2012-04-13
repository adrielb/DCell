#include "petsc.h"
#include "petscda.h"
#include "Common.h"
#include <unistd.h>
#include <mpi.h>

PetscErrorCode AdvectSL( DA daVel, Vec vecVel, DA daBuf, Vec vecBuf, DA daPhi, Vec vecPhi, PetscReal *phiShift, int f, PetscReal dh[3], PetscReal dt );
PetscErrorCode AdvectSL_2D( DA daVel, Vec vecVel, DA daBuf, Vec vecBuf, DA daPhi, Vec vecPhi, PetscReal *phiShift, int f, PetscReal dh[3], PetscReal dt );

const char* prefix = "/home/abergman/Research/DCell/temp/";
//const char* prefix = "/scratch/temp/";

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc, char **args)
{
  PetscLogDouble t1,t2;
  PetscErrorCode  ierr;
  ierr = PetscInitialize(&argc, &args, (char *) 0, ""); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "Start %s\n", __FILE__); CHKERRQ(ierr);
  MPI_Comm comm = PETSC_COMM_WORLD;

  int rank;
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
  char hostname[256];
  gethostname(hostname,255);
  printf("[%d] %s\n",rank,hostname);
  PetscBarrier(0);

  int MX = 64;
  int dof = 2;

  DA da;
  ierr = DACreate2d(PETSC_COMM_WORLD,DA_NONPERIODIC,DA_STENCIL_BOX,
          MX,MX,  PETSC_DECIDE,PETSC_DECIDE,  dof,1,  0,0, &da); CHKERRQ(ierr);
  ierr = DAView(da,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);

  Vec vel;
  ierr = DACreateGlobalVector(da,&vel); CHKERRQ(ierr);
  PetscReal ***v;
  ierr = DAVecGetArrayDOF(da,vel,&v); CHKERRQ(ierr);
  int i,j;
  DALocalInfo a;
  ierr = DAGetLocalInfo(da,&a); CHKERRQ(ierr);
  PetscReal dx=1/(MX-2.),dy=1/(MX-2.),dz=1/(MX-2.);
  PetscReal maxVel = 0, curVel = 0;
  for (j = a.ys; j < a.ys+a.ym; ++j) {
    for (i = a.xs; i < a.xs+a.xm; ++i) {
      v[j][i][0] = -(j+0.5 - a.mx/2.) * dx;
      v[j][i][1] = -(i+0.5 - a.my/2.) * dy;
//      v[j][i][1] = 1;
//      v[j][i][1] = 1;

      curVel = sqrt( v[j][i][0]*v[j][i][0] +
                     v[j][i][1]*v[j][i][1] );
      maxVel = maxVel < curVel ? curVel : maxVel;
    }
  }
  ierr = DAVecRestoreArray(da,vel,&v); CHKERRQ(ierr);

  ierr = MPI_Allreduce(&maxVel,&maxVel,1,MPI_DOUBLE,MPI_MAX,comm); CHKERRQ(ierr);

  char filename[128];
  PetscViewer binv;
  sprintf(filename,"%s/uvp.Real64",prefix);
  ierr = PetscViewerBinaryOpen(comm,filename,FILE_MODE_WRITE,&binv); CHKERRQ(ierr);
  ierr = VecView(vel, binv); CHKERRQ(ierr);
  ierr = PetscViewerDestroy(binv); CHKERRQ(ierr);

  DA daPhi, daBuf;
  ierr = DACreate2d(comm,DA_NONPERIODIC,DA_STENCIL_BOX, MX,MX, PETSC_DECIDE,PETSC_DECIDE, 2,1, 0,0, &daPhi); CHKERRQ(ierr);
  ierr = DACreate2d(comm,DA_NONPERIODIC,DA_STENCIL_BOX, MX,MX, PETSC_DECIDE,PETSC_DECIDE, 1,1, 0,0, &daBuf); CHKERRQ(ierr);

  Vec vecPhi, vecBuf;
  ierr = DACreateGlobalVector(daPhi,&vecPhi); CHKERRQ(ierr);
  ierr = DACreateGlobalVector(daBuf,&vecBuf); CHKERRQ(ierr);

  // phi
  PetscReal ***phi;
  ierr = DAVecGetArrayDOF(daPhi,vecPhi,&phi); CHKERRQ(ierr);
  for   (j = a.ys; j < a.ys+a.ym; ++j) {
    for (i = a.xs; i < a.xs+a.xm; ++i) {
      if ( 0.15 < i*dx && i*dx < 0.85 &&
           0.45 < j*dy && j*dy < 0.55  )
      {
        phi[j][i][0] = 1;
      }
    }
  }
  ierr = DAVecRestoreArray(daPhi,vecPhi,&phi); CHKERRQ(ierr);

  PetscReal CFL = 1;
  PetscReal dt = CFL * dx / maxVel;
  PetscReal dh[3] = {dx,dy,dz};
  PetscReal phiShift[3] = {0,0,0};

  int t;
  PetscPrintf(comm,"maxVel: %f\n", maxVel);
  PetscPrintf(comm,"dh: (%f,%f,%f)\n", dh[0],dh[1],dh[2]);
  PetscPrintf(comm,"dt: %f\n", dt);

  sprintf(filename,"%s/phi.100.Real64",prefix);
  ierr = PetscViewerBinaryOpen(comm,filename,FILE_MODE_WRITE,&binv); CHKERRQ(ierr);
  ierr = VecView(vecPhi, binv); CHKERRQ(ierr);
  ierr = PetscViewerDestroy(binv); CHKERRQ(ierr);

  for (t = 101; t < 205; ++t) {

    ierr = AdvectSL_2D(da,vel,daBuf,vecBuf,daPhi,vecPhi,phiShift,0,dh,dt); CHKERRQ(ierr);

    if( t == 150 )
    {
      ierr = VecScale(vel,-1); CHKERRQ(ierr);
    }

    ierr = PetscGetTime(&t1); CHKERRQ(ierr);
    sprintf(filename,"%s/phi.%d.Real64",prefix,t);
    ierr = PetscViewerBinaryOpen(comm,filename,FILE_MODE_WRITE,&binv); CHKERRQ(ierr);
    ierr = VecView(vecPhi, binv); CHKERRQ(ierr);
    ierr = PetscViewerDestroy(binv); CHKERRQ(ierr);
    ierr = PetscGetTime(&t2); CHKERRQ(ierr);
    ierr = PetscPrintf(comm,"Vec View: %f sec\n",t2-t1); CHKERRQ(ierr);
  }

  ierr = PetscPrintf(PETSC_COMM_WORLD, "End %s\n", __FILE__); CHKERRQ(ierr);
  ierr = PetscFinalize(); CHKERRQ(ierr);
}

PetscErrorCode test3D() {
  PetscErrorCode  ierr;
  MPI_Comm comm = PETSC_COMM_WORLD;
  int MX = 64;
  int dof = 2;

  DA da;
  ierr = DACreate3d(PETSC_COMM_WORLD,DA_NONPERIODIC,DA_STENCIL_STAR,
          MX,MX,MX,
          PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE,
          dof,1,
          0,0,0, &da); CHKERRQ(ierr);
  ierr = DAView(da,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);

  Vec vel;
  ierr = DACreateGlobalVector(da,&vel); CHKERRQ(ierr);
  PetscReal ****v;
  ierr = DAVecGetArrayDOF(da,vel,&v); CHKERRQ(ierr);
  int i,j,k;
  int xs,ys,zs,xm,ym,zm;
  ierr = DAGetCorners(da,&xs,&ys,&zs,&xm,&ym,&zm); CHKERRQ(ierr);
  PetscReal p = PETSC_PI, p2 = 2*PETSC_PI;
  PetscReal Ux,Uy,Uz, Vx,Vy,Vz, Wx,Wy,Wz;
  PetscReal dx=1/(MX-2.),dy=1/(MX-2.),dz=1/(MX-2.);
  PetscReal maxVel = 0, curVel = 0;
  for (k = zs; k < zm; ++k) {
    for (j = ys; j < ym; ++j) {
      for (i = xs; i < xm; ++i) {
        Ux = (i-0.5)*dx; Uy = j*dy;       Uz = dz*k;
        Vx = i*dx;       Vy = (j-0.5)*dy; Vz = dz*k;
        Wx = i*dx;       Wy = j*dy;       Wz = dz*(k-0.5);
        v[k][j][i][0] = 2*sin(p*Ux)*sin(p*Ux) * sin(p2*Uy) * sin(p2*Uz);
        v[k][j][i][1] =  -sin(p*Vy)*sin(p*Vy) * sin(p2*Vx) * sin(p2*Vz);
        v[k][j][i][2] =  -sin(p*Wz)*sin(p*Wz) * sin(p2*Wx) * sin(p2*Wy);
        curVel = sqrt( v[k][j][i][0]*v[k][j][i][0] +
                       v[k][j][i][1]*v[k][j][i][1] +
                       v[k][j][i][2]*v[k][j][i][2] );
        maxVel = maxVel < curVel ? curVel : maxVel;
      }
    }
//    printf("k: %f\n", Wz);
  }
  /* Divergence check
  for (k = zs; k < zm-1; ++k) {
    for (j = ys; j < ym-1; ++j) {
      for (i = xs; i < xm-1; ++i) {
        printf("%f  ", v[k][j][i+1][0] - v[k][j][i][0] +
                       v[k][j+1][i][1] - v[k][j][i][1] +
                       v[k+1][j][i][2] - v[k][j][i][2] );
      }
      printf("\n");
    }
  }*/
  ierr = DAVecRestoreArray(da,vel,&v); CHKERRQ(ierr);

  PetscViewer binv;
  ierr = PetscViewerBinaryOpen(comm,"/scratch/temp/uvp.Real64",FILE_MODE_WRITE,&binv); CHKERRQ(ierr);
  ierr = VecView(vel, binv); CHKERRQ(ierr);
  ierr = PetscViewerDestroy(binv); CHKERRQ(ierr);


  DA daPhi, daBuf;
  ierr = DACreate3d(comm,DA_NONPERIODIC,DA_STENCIL_STAR, MX,MX,MX, PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE, 2,1, 0,0,0, &daPhi); CHKERRQ(ierr);
  ierr = DACreate3d(comm,DA_NONPERIODIC,DA_STENCIL_STAR, MX,MX,MX, PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE, 1,1, 0,0,0, &daBuf); CHKERRQ(ierr);

  Vec vecPhi, vecBuf;
  ierr = DACreateGlobalVector(daPhi,&vecPhi); CHKERRQ(ierr);
  ierr = DACreateGlobalVector(daBuf,&vecBuf); CHKERRQ(ierr);

  // phi
  PetscReal ****phi;
  ierr = DAVecGetArrayDOF(daPhi,vecPhi,&phi); CHKERRQ(ierr);
  for (k = zs; k < zm; ++k) {
    for (j = ys; j < ym; ++j) {
      for (i = xs; i < xm; ++i) {
        if ( 20 < i && i < 30 &&
             20 < j && j < 30 &&
             20 < k && k < 30  )
        {
          phi[k][j][i][0] = 1;
        }
      }
    }
  }
/*
  PetscReal shift[3] = {0,0,0}, X[3], a;
  for (k = zs; k < zm-1; ++k) {
    for (j = ys; j < ym-1; ++j) {
      for (i = xs; i < xm-1; ++i) {
        X[0] = i+0.5; X[1] = j+0.5; X[2] = k+0.5;
//        X[0] = i; X[1] = j; X[2] = k+0.5;
        a = InterpolateField3D( phi, 0, shift, X);
        printf("%1.1f ",a);
      }
      printf("\n");
    }
    printf("\n");
  }
  exit(0);
*/
  ierr = DAVecRestoreArray(daPhi,vecPhi,&phi); CHKERRQ(ierr);

}
