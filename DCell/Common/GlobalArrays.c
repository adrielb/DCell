#include "Common.h"

PetscLogEvent EVENT_GAPutVec;
PetscLogEvent EVENT_GAGetVec;
PetscLogEvent EVENT_GAGather;
PetscLogEvent EVENT_GAScatterAcc;


#undef __FUNCT__
#define __FUNCT__ "GACreate"
PetscErrorCode GACreate( DM da, int *ga )
{
  int i;
  int dims[4], block[4];
  int dim, M,N,P, m,n,p, dof;
  const PetscInt *lx, *ly, *lz;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = DMDAGetOwnershipRanges(da,&lx,&ly,&lz); CHKERRQ(ierr);
  ierr = DMDAGetInfo(da,&dim,&M,&N,&P,&m,&n,&p,&dof,0,0,0,0,0); CHKERRQ(ierr);
  dim = dim + 1; // plus one for [dof]
  p = p < 0 ? 0 : p;
  int map[p+n+m+1]; // size of map is sum of block[]
  
  if( dim == 3 ) {  // 2D field
    dims[0] = N;   block[0] = n;
    dims[1] = M;   block[1] = m;
    dims[2] = dof; block[2] = 1;

    map[0] = 0;
    for( i = 1; i < n; i++ ) map[i] = ly[i-1] + map[i-1];
    map[n] = 0;
    for( i = 1; i < m; i++ ) map[i+n] = lx[i-1] + map[i+n-1];
    map[n+m] = 0;

  } else { // 3D field
    dims[0] = P;   block[0] = p;
    dims[1] = N;   block[1] = n;
    dims[2] = M;   block[2] = m;
    dims[3] = dof; block[3] = 1;
    
    map[0] = 0;
    for( i = 1; i < p; i++ ) map[i] = lz[i-1] + map[i-1];
    map[p] = 0;
    for( i = 1; i < n; i++ ) map[i+p] = ly[i-1] + map[i+p-1];
    map[p+n] = 0;
    for( i = 1; i < m; i++ ) map[i+p+n] = lx[i-1] + map[i+p+n-1];
    map[p+n+m] = 0;
  }

  int g = GA_Create_handle();
  GA_Set_data(g,dim,dims,C_DBL);
  GA_Set_irreg_distr(g,map,block);
  ierr = GA_Allocate( g ); if( !ierr ) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_MEM,"GA_Allocate() Failed");
  *ga = g;
  
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "GAPutVec"
PetscErrorCode GAPutVec(Vec vec, int ga )
{
  int len;
  int type,dim,dims[4];
  int nodeid = GA_Nodeid();
  PetscReal *petscvec;
  PetscErrorCode ierr;
  PetscFunctionBegin;
  ierr = PetscLogEventBegin(EVENT_GAPutVec,0,0,0,0); CHKERRQ(ierr);
  NGA_Inquire(ga,&type,&dim,dims);
  int lo[dim],hi[dim],ld[3];
  ierr = VecGetSize(vec,&len); CHKERRQ(ierr);
  ierr = VecGetArray(vec,&petscvec); CHKERRQ(ierr);
  NGA_Distribution(ga,nodeid,lo,hi);
  ld[0] = hi[1]-lo[1]+1;
  ld[1] = hi[2]-lo[2]+1;
  ld[2] = dim==4 ? hi[3]-lo[3]+1 : 0;
  NGA_Put(ga,lo,hi,petscvec,ld);
  ierr = VecRestoreArray(vec,&petscvec); CHKERRQ(ierr);
  GA_Sync();
  ierr = PetscLogEventEnd(EVENT_GAPutVec,0,0,0,0); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "GAGetVec"
PetscErrorCode GAGetVec(int ga, Vec vec )
{
  int len;
  int type,dim,dims[4];
  int nodeid = GA_Nodeid();
  PetscReal *petscvec;
  PetscErrorCode ierr;
  PetscFunctionBegin;
  ierr = PetscLogEventBegin(EVENT_GAGetVec,0,0,0,0); CHKERRQ(ierr);
  NGA_Inquire(ga,&type,&dim,dims);
  int lo[dim],hi[dim],ld[3];
  ierr = VecGetSize(vec,&len); CHKERRQ(ierr);
  ierr = VecGetArray(vec,&petscvec); CHKERRQ(ierr);
  NGA_Distribution(ga,nodeid,lo,hi);
  ld[0] = hi[1]-lo[1]+1;
  ld[1] = hi[2]-lo[2]+1;
  ld[2] = dim==4 ? hi[3]-lo[3]+1 : 0;
  NGA_Get(ga,lo,hi,petscvec,ld);
  ierr = VecRestoreArray(vec,&petscvec); CHKERRQ(ierr);
  GA_Sync();
  ierr = PetscLogEventEnd(EVENT_GAGetVec,0,0,0,0); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "GAGather"
PetscErrorCode GAGather(int ga, Array c, Array v)
{
  PetscReal *val;
  int **coor = ArrayGetData(c);
  int len = ArrayLength(c);
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscLogEventBegin(EVENT_GAGather,0,0,0,0); CHKERRQ(ierr);
  ierr = ArraySetSize(v,len); CHKERRQ(ierr);
  val =  ArrayGetData(v);
  NGA_Gather(ga,val,coor,len);
  ierr = PetscLogEventEnd(EVENT_GAGather,0,0,0,0); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "GAScatterAcc"
PetscErrorCode GAScatterAcc(int ga, Array c, Array v)
{
  PetscReal alpha = 1;
  PetscReal *val = ArrayGetData(v);
  int **coor = ArrayGetData(c);
  int len = ArrayLength(c);
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscLogEventBegin(EVENT_GAScatterAcc,0,0,0,0); CHKERRQ(ierr);
  NGA_Scatter_acc(ga,val,coor,len,&alpha);
  ierr = ArraySetSize(c,0); CHKERRQ(ierr);
  ierr = ArraySetSize(v,0); CHKERRQ(ierr);
  ierr = PetscLogEventEnd(EVENT_GAScatterAcc,0,0,0,0); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
