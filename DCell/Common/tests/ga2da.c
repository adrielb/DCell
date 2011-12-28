#include "Common.h"

const PetscInt dof = 3;
DA da;
DALocalInfo info;
int rank;

#undef __FUNCT__
#define __FUNCT__ "TestVecs"
PetscErrorCode TestVecs( int ga, Vec petscVec)
{
  PetscErrorCode  ierr;
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Testing GA_Get()\t"); CHKERRQ(ierr);
  PetscReal ***petscvec;
  int c = 0;
  int lo[3],hi[3],ld[2]={info.xm,dof};
  int nodeid = GA_Nodeid();
  PetscReal *gabuf;
  ierr = PetscMalloc(sizeof(PetscReal)*info.xm*info.ym*dof,&gabuf); CHKERRQ(ierr);
  NGA_Distribution(ga,nodeid,lo,hi);
  NGA_Get(ga,lo,hi,gabuf,ld);
//  ierr = PetscSynchronizedPrintf(PETSC_COMM_WORLD,"GA[%d]: lo [%d,%d,%d]   hi [%d,%d,%d]   ld [%d,%d]\n",nodeid, lo[0],lo[1],lo[2],hi[0],hi[1],hi[2],ld[0],ld[1]); CHKERRQ(ierr);
//  ierr = PetscSynchronizedFlush(PETSC_COMM_WORLD); CHKERRQ(ierr);
  ierr = DAVecGetArrayDOF(da,petscVec,&petscvec); CHKERRQ(ierr);
  int j,i,d;
  for ( j = info.ys; j < info.ys+info.ym; ++j) {
    for ( i = info.xs; i < info.xs+info.xm; ++i) {
//      printf("%d:(%d,%d)", rank, i, j);
      for ( d = 0; d < dof; ++d) {
//        printf("\t\t[%1.0f,%1.0f]", petscvec[j][i][d], gabuf[c] );
        if( petscvec[j][i][d] != gabuf[c] ) {
          printf("[%d,%d;%d] petsc: %1.0f, ga: %1.0f\n",j,i,d,petscvec[j][i][d],gabuf[c]);
          SETERRQ1(0,"[%d]GA_Get() FAILED", rank );
        }
        c++;
      }
//      printf("\n");
    }
//    printf("\n");
  }
  ierr = DAVecRestoreArrayDOF(da,petscVec,&petscvec); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"PASSED [%d]\n",rank); CHKERRQ(ierr);
  PetscBarrier(0);


  ierr = PetscPrintf(PETSC_COMM_WORLD,"Testing GA_Gather()\n"); CHKERRQ(ierr);
  if( rank == 0 ) {
    printf("Testing Gather_flat (single element)...");
    PetscReal value;
    int subscript[dof];
    c=0;
    for ( j = 0; j < info.my; ++j) {
      for ( i = 0; i < info.mx; ++i) {
        for ( d = 0; d < dof; ++d) {
          subscript[2] = d;
          subscript[1] = i;
          subscript[0] = j;
          NGA_Gather_flat(ga,&value,subscript,1);
          if( PetscAbs( value - c/3 ) > 1e-3 ) {
            printf( "[ %d, %d; %d] ", i, j, d);
            printf("%1.1f - %d \n", value, c/3);
            SETERRQ(0,"GA_Gather() FAILED");
          }
          c++;
        }
      }
    }
    printf("Success\n");
    printf("Testing Gather_flat (dof element)...");
    int subscriptn[3*dof];
    PetscReal values[3];
    c=0;
    for ( j = 0; j < info.my; ++j) {
      for ( i = 0; i < info.mx; ++i) {
        for ( d = 0; d < dof; ++d) {
          subscriptn[d*dof+2] = d;
          subscriptn[d*dof+1] = i;
          subscriptn[d*dof+0] = j;
        }
        NGA_Gather_flat(ga,values,subscriptn,dof);
        printf( "[ %d, %d] ", i, j);
        printf("%1.1f %1.1f %1.1f %d \n", values[0], values[1], values[2], c);
        c++;
      }
    }
    printf("Success\n");
    
    int n = info.mx*info.my*dof;
    PetscReal *val;
    int *coor, **coor2;
    PetscMalloc(sizeof(PetscReal)*n,&val);
    PetscMalloc(sizeof(int) *n*3,&coor);
    PetscMalloc(sizeof(int*)*n,&coor2);
    for ( i = 0; i < n; ++i) {
      coor2[i] = &coor[3*i];
    }
    c=0;
    for ( j = 0; j < info.my; ++j) {
      for ( i = 0; i < info.mx; ++i) {
        for ( d = 0; d < dof; ++d) {
          coor2[c][0] = d;
          coor2[c][1] = i;
          coor2[c][2] = j;
          c++;
        }
      }
    }
    /* printing coor
    for (int i = 0; i < n;i++) {
      printf("[%d,%d;%d] ",coor2[i][0],coor2[i][1],coor2[i][2]);
      if( i%3 == 2 ) printf("\n");
    }
    printf("\n\n\n");
    */
	printf("\n");
    NGA_Gather(ga,val,coor2,n);
    for ( i = 0; i < n; ++i) {
//      if( PetscAbs(val[i] - i/3) > 1e-3 ) {
//        printf("[%d,%d;%d] petsc: %d, ga: %1.0f\n",coor2[i][2],coor2[i][1],coor2[i][0],i/3,val[i]);
//        SETERRQ(0,"GA_Gather() FAILED");
//      }
      printf("(%d,%d;%d): %1.0f \t",coor2[i][0],coor2[i][1],coor2[i][2],val[i]);
      if( i%3 == 2 ) printf("\n");
    }
    NGA_Scatter(ga,val,coor2,n);
  }
  PetscBarrier(0);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"PASSED\n"); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode MakeGAVec( int d1, int d2, int *p_ga ) {
  PetscErrorCode ierr;
  int dim, M,N,P, m,n,p, dof;
  ierr = DAGetInfo(da,&dim,&M,&N,&P,&m,&n,&p,&dof,0,0,0); CHKERRQ(ierr);
  dim = dim + 1; // plus one for [dof]
  int dims[3] = {N,M,dof};
  int block[3] = {n,m,1};
  int map[n+m+1]; // size of map is sum of block[]

  const PetscInt *lx, *ly, *lz;
  DAGetOwnershipRanges(da,&lx,&ly,&lz);
  map[0] = 0;
  int i;
  for( i = 1; i < n; i++ ) map[i] = ly[i-1] + map[i-1];
  map[n] = 0;
//  for( i = n+1; i < m+n; i++ ) map[i] = lx[i-n-1] + map[i-1];
  for( i = 1; i < m; i++ ) map[i+n] = lx[i-1] + map[i+n-1];
  map[n+m] = 0;

  if( rank == 0 ) {
    for ( i = 0; i < m; ++i) {
      printf("lx[%d]= %d\n", i, lx[i]);
    }
    for ( i = 0; i < n; ++i) {
      printf("ly[%d]= %d\n", i, ly[i]);
    }
    for ( i = 0; i < p; ++i) {
      printf("lz[%d]= %d\n", i, lz[i]);
    }
    for (i = 0; i < n+m+1; ++i) {
      printf("map[%d]=%d\n",i,map[i]);
    }
    for ( i = 0; i < dim; ++i) {
      printf("block[%d]=%d \t dims[%d]=%d\n",i,block[i],i,dims[i]);
    }
  }
  int ga = GA_Create_handle();
  GA_Set_data(ga,dim,dims,C_DBL);
  GA_Set_irreg_distr(ga,map,block);
  ierr = GA_Allocate( ga ); if( !ierr ) exit(1);
//  int ga = NGA_Create_irreg(C_DBL,dim,dims,name,map,block);
//  int ga = NGA_Create(C_DBL,dim,dims,name,0);
  int nodeid = GA_Nodeid();
  printf("GA: %d   MPI: %d\n", nodeid, rank);
  int low[dim], hi[dim];
  NGA_Distribution(ga,nodeid,low,hi);
  ierr = PetscSynchronizedPrintf(PETSC_COMM_WORLD,"GA[%d]: low [%d,%d,%d]   hi [%d,%d,%d]\n",nodeid, low[0],low[1],low[2],hi[0],hi[1],hi[2]); CHKERRQ(ierr);
  ierr = PetscSynchronizedFlush(PETSC_COMM_WORLD); CHKERRQ(ierr);
//  GA_Print(ga);
  *p_ga = ga;
  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "MakePetscVec"
PetscErrorCode MakePetscVec( int d1, int d2, Vec *vec)
{
  const PetscInt sw = 1;
  const PetscInt m = PETSC_DECIDE, n = PETSC_DECIDE;
  Vec petscVec;
  PetscReal ***petscvec;

  PetscErrorCode  ierr;
  ierr = DACreate2d(PETSC_COMM_WORLD,DA_NONPERIODIC,DA_STENCIL_STAR,
           d1,d2,m,n,dof,sw,0,0, &da); CHKERRQ(ierr);
  ierr = DAGetLocalInfo(da,&info); CHKERRQ(ierr);
  ierr = DAView(da, PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
  ierr = DACreateGlobalVector(da, &petscVec); CHKERRQ(ierr);
  ierr = DAVecGetArrayDOF(da,petscVec,&petscvec); CHKERRQ(ierr);
  int j,i,d;
  for ( j = info.ys; j < info.ys+info.ym; ++j) {
    for ( i = info.xs; i < info.xs+info.xm; ++i) {
      for ( d = 0; d < dof; ++d) {
        petscvec[j][i][d] = i + j * d1;
      }
    }
  }
  ierr = DAVecRestoreArrayDOF(da,petscVec,&petscvec); CHKERRQ(ierr);
  ierr = PetscBarrier(0); CHKERRQ(ierr);
//  ierr = VecView(petscVec,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);

  *vec = petscVec;
  PetscFunctionReturn(0);
}

int main(int argc, char **args)
{
  PetscErrorCode  ierr;
  ierr = PetscInitialize(&argc, &args, (char *) 0, ""); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "Start %s\n", __FILE__); CHKERRQ(ierr);

  int stack = 100000000;
  int heap =  10000000;
  if( !MA_init(C_DBL, stack, heap)) {
    GA_Error((char*)"MA_init failed",stack+heap);
  }
  GA_Initialize();

  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
  const size_t len = 1024;
  char name[len];
  PetscGetHostName(name,len);
  ierr = PetscSynchronizedPrintf(PETSC_COMM_WORLD,"%d: %s\n", rank, name); CHKERRQ(ierr);
  ierr = PetscSynchronizedFlush(PETSC_COMM_WORLD); CHKERRQ(ierr);

  int d1 = 5, d2 = d1;

  int gavec;
  Vec petscvec;
  ierr = MakePetscVec(d1,d2,&petscvec); CHKERRQ(ierr);
  ierr = GACreate(da,&gavec); CHKERRQ(ierr);
  ierr = GAPutVec(petscvec,gavec); CHKERRQ(ierr);
  ierr = TestVecs(gavec,petscvec); CHKERRQ(ierr);
  ierr = GAGetVec(gavec,petscvec); CHKERRQ(ierr);
  GA_Print_distribution(gavec);
  GA_Destroy(gavec);
  if(rank==0) {
    GA_Print_stats();
    GA_Summarize(1);
  }


//  GA_Uses_ma() ? printf("GA uses MA\n") : printf("Not using MA\n");

  GA_Terminate();
  ierr = PetscPrintf(PETSC_COMM_WORLD, "End %s\n", __FILE__); CHKERRQ(ierr);
  ierr = PetscFinalize(); CHKERRQ(ierr);
}
