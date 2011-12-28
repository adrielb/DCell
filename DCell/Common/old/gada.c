#include "petsc.h"
#include "petscda.h"
#include "ga.h"
#include "macdecls.h"
//#include "/share/apps/petsc-3.1-p3/include/private/daimpl.h" //TODO: change include file to relative path

PetscErrorCode DACreateGlobalArray( DA da, int *GA, Vec *g );

#undef __FUNCT__
#define __FUNCT__ "testGet"
PetscErrorCode testGet(  )
{
  PetscErrorCode ierr;
  
  
  PetscFunctionBegin;
  
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "testCreate3D"
PetscErrorCode testCreate3D(  )
{
  int ga;
  DA da;
  DALocalInfo info;
  Vec vec;
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  int d1 = 229, d2 = 229, d3 = 229;
  int rank;
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
  ierr = DACreate3d(PETSC_COMM_WORLD,DA_NONPERIODIC,DA_STENCIL_STAR,
              d1,d2,d3,
              PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE,
              1,1,
              0,0,0, &da); CHKERRQ(ierr);
  ierr = DAGetLocalInfo(da,&info); CHKERRQ(ierr);
//  ierr = DACreateGlobalArray( da, &ga, &vec); CHKERRQ(ierr);
  
  PetscReal ***v;
  ierr = DAVecGetArray(da,vec,&v); CHKERRQ(ierr);
  int xe = info.xs+info.xm,
      ye = info.ys+info.ym,
      ze = info.zs+info.zm;
  int i,j,k;
  for (k = info.zs; k < ze; ++k) {
    for (j = info.ys; j < ye; ++j) {
      for (i = info.xs; i < xe; ++i) {
        v[k][j][i] = 1.*i + d1*j + d1*d2*k;
      }
    }
  }
  ierr = DAVecRestoreArray(da,vec,&v); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "Sequential values filled in petsc vec.\n"); CHKERRQ(ierr);
  
  ierr = PetscBarrier(0); CHKERRQ(ierr);
  int lo[3],ld, p = 10;
  int patch[10][10][10];
  double val;
  for (int k = 0; k < d3; k+=p) {
    for (int j = 0; j < d2; j+=p) {
      for (int i = 0; i < d1; i+=p) {
        lo[0] = k;
        lo[1] = j;
        lo[2] = i;
        NGA_Get(ga,lo,lo,&val,&ld);
        if( PetscAbs( i + d1*j + d1*d2*k - val) > .1 )
//          printf(".");
          printf("(%3.0f,%3.0f) ", 1.*i + d1*j + d1*d2*k, val);
      }
    }
  }
  ierr = PetscPrintf(PETSC_COMM_WORLD, "Ended NGA_Get() test.\n"); CHKERRQ(ierr);
  
  ierr = PetscBarrier(0); CHKERRQ(ierr);
  if( rank == 0 )
  {
    for (int k = 0; k < d3; ++k) {
      printf(">%d\n",k);
      for (int j = 0; j < d2; ++j) {
        for (int i = 0; i < d1; ++i) {
          lo[0] = k; lo[1] = j; lo[2] = i;
          val = 1.*i + d1*j + d1*d2*k;
          val *= -1;
          NGA_Put(ga,lo,lo,&val,&ld);
        }
      }
    }
  }
  ierr = PetscPrintf(PETSC_COMM_WORLD, "Ended NGA_Put() negative seq values.\n"); CHKERRQ(ierr);
  
  ierr = PetscBarrier(0); CHKERRQ(ierr);
  ierr = DAVecGetArray(da,vec,&v); CHKERRQ(ierr);
  for (int k = info.zs; k < ze; ++k) {
    for (int j = info.ys; j < ye; ++j) {
      for (int i = info.xs; i < xe; ++i) {
        val = -1 * (1.*i + d1*j + d1*d2*k);
        if( PetscAbs( val - v[k][j][i] ) > .1 )
          printf(".");
      }
    }
  }
  ierr = DAVecRestoreArray(da,vec,&v); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "Ended petsc vec update test.\n"); CHKERRQ(ierr);
  
  if( rank == 0 )
    GA_Print_stats();
  
  ierr = VecDestroy(vec); CHKERRQ(ierr);
  GA_Destroy(ga);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "testCreate2D"
PetscErrorCode testCreate2D()
{
  int ga;
  DA da;
  DALocalInfo info;
  Vec vec;
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  int d1 = 1453, d2 = 1451;
  ierr = DACreate2d(PETSC_COMM_WORLD,DA_NONPERIODIC,DA_STENCIL_STAR,
              d1,d2,PETSC_DECIDE,PETSC_DECIDE,1,1,0,0, &da); CHKERRQ(ierr);
  ierr = DAGetLocalInfo(da,&info); CHKERRQ(ierr);
//  ierr = DACreateGlobalArray( da, &ga, &vec); CHKERRQ(ierr);
  
  PetscReal **v;
  ierr = DAVecGetArray(da,vec,&v); CHKERRQ(ierr);
  int xe = info.xs+info.xm,
      ye = info.ys+info.ym;
  for (int j = info.ys; j < ye; ++j) {
    for (int i = info.xs; i < xe; ++i) {
      v[j][i] = 1.*i + d1 * j;
    }
  }
  ierr = DAVecRestoreArray(da,vec,&v); CHKERRQ(ierr);
  PetscPrintf(PETSC_COMM_WORLD,"Updated local portion with DAVec\n");
  PetscBarrier(0);
  {
    double *da_ptr;
  VecGetArray(vec, &da_ptr);
  double *ptr;
  int low[2],hi[2],ld;
  NGA_Distribution(ga,GA_Nodeid(),low,hi);
  NGA_Access(ga,low,hi,&ptr,&ld);
  printf("[%d] ga:%p\tda:%p\tdiff:%p\n", GA_Nodeid(), ptr, da_ptr, (ptr-da_ptr) );
  NGA_Release_update(ga,low,hi);
  }
  
  int lo[2],ld;
  double val;
  for (int j = 0; j < d2; ++j) {
    for (int i = 0; i < d1; ++i) {
      lo[0] = j;
      lo[1] = i;
      NGA_Get(ga,lo,lo,&val,&ld);
      if( PetscAbs( i + d1*j - val) > .1 )
        printf(".");
//        printf("[%d] (%3.0f,%3.0f)\n", GA_Nodeid(), 1.*i + d1*j, val);
    }
  }
  GA_Print_stats();
  ierr = VecDestroy(vec); CHKERRQ(ierr);
  GA_Destroy(ga);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "vizGA2DA"
PetscErrorCode vizGA2DA()
{
  PetscErrorCode  ierr;
  int rank;
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);  
  int d1 = 40, d2 = 50;
  
  DA da;
  Vec vec;
  const PetscInt *lx, *ly, *lz;
  const PetscInt dof = 1, sw = 1;
  PetscInt m=1,n=1,p;
  DALocalInfo info;
  ierr = DACreate2d(PETSC_COMM_WORLD,DA_NONPERIODIC,DA_STENCIL_STAR,
         d1,d2,m,n,dof,sw,0,0, &da); CHKERRQ(ierr);
  ierr = DACreateGlobalVector(da, &vec); CHKERRQ(ierr);
  ierr = DAGetOwnershipRanges(da, &lx, &ly, &lz); CHKERRQ(ierr);
  ierr = DAGetLocalInfo(da,&info); CHKERRQ(ierr);
  ierr = DAGetInfo(da,0,0,0,0,&m,&n,&p,0,0,0,0); CHKERRQ(ierr);
  /**/
  ierr = DAView(da, PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
  for (int i = 0; i < m; ++i) {
    PetscPrintf(PETSC_COMM_WORLD,"%d\tlx: %d\n",i,lx[i]);
  }
  for (int i = 0; i < n; ++i) {
    PetscPrintf(PETSC_COMM_WORLD,"%d\tly: %d\n",i,ly[i]);
  }
  /**/
 
  exit(0);
  int ga = GA_Create_handle();
  int ndim = 2;
  int dims[2] = {d2,d1};
  GA_Set_data(ga,2,dims,MT_F_DBL);
  int *map;
  PetscMalloc( sizeof(int)*(m+n), &map);
  map[0] = 0;
  for( int i = 1; i < n; i++ )
  {
    map[i] = ly[i-1] + map[i-1];
  }
  map[n] = 0;
  for( int i = n+1; i < m+n; i++ )
  {
    map[i] = lx[i-n-1] + map[i-1];
  }
  /* correct ordering, but nodeid's dont line up with mpi rank for petsc's da
   * DA: +---+---+   GA: +---+---+   
   *     +-2-+-3-+       +-1-+-3-+
   *     +---+---+       +---+---+
   *     +-0-+-1-+       +-0-+-2-+
   *     +---+---+       +---+---+
  int *map;
  PetscMalloc( sizeof(int)*(m+n), &map);
  map[0] = 0;
  for( int i = 1; i < m; i++ )
  {
    map[i] = lx[i] + map[i-1];
  }
  map[m] = 0;
  for( int i = m+1; i < m+n; i++ )
  {
    map[i] = ly[i-m] + map[i-1];
  }
  */
  int block[2] = {n,m};  
  GA_Set_irreg_distr(ga,map,block);
  ierr = GA_Allocate( ga );
  if( !ierr ) GA_Error("\n\n\nga allocaltion failed\n\n",ierr);
  if( !ga ) GA_Error("\n\n\n ga null \n\n",ierr); 
  if( rank != GA_Nodeid() ) GA_Error("MPI rank does not match GA_Nodeid()",1);
  GA_Print_distribution(ga);  
  
  int lo[2], hi[2];
  NGA_Distribution(ga,rank,lo,hi);
  if( lo[1] != info.xs || hi[1] != info.xs+info.xm-1 ||
      lo[0] != info.ys || hi[0] != info.ys+info.ym-1 )
  {
    PetscSynchronizedPrintf(PETSC_COMM_SELF,"[%d] lo:(%2d,%2d)  hi:(%2d,%2d) \t DA: (%2d,%2d), (%2d, %2d)\n",
        rank, lo[1], lo[0], hi[1], hi[0], info.xs, info.ys, info.xs+info.xm-1, info.ys+info.ym-1);
  }
  PetscBarrier(0);
  PetscSynchronizedFlush(PETSC_COMM_WORLD);

  AO ao;
  DAGetAO(da,&ao);
  if( rank == 0 )
  {
    int *idx, len = d1*d2;
    PetscReal *val;
    PetscMalloc(sizeof(PetscReal)*len, &val);
    PetscMalloc(sizeof(int)*len, &idx);
    for (int j = 0; j < d2; ++j)
    {
      for (int i = 0; i < d1; ++i)
      {
        idx[i + d1*j] = i + d1*j;
        val[i + d1*j] = i + d1*j;
      }
    }
    AOApplicationToPetsc(ao,len,idx);
    VecSetValues(vec,len,idx,val,INSERT_VALUES);

    int a[2], b[2],ld[1]={0};
    double c = 0;
    for (int j = 0; j < d2; ++j)
    {
      for (int i = 0; i < d1; ++i)
      {
        a[0] = j;
        a[1] = i;
//        printf("%5.0f ",c);
        NGA_Put(ga,a,a,&c,ld);
        c++;
      }
    }
  }
//  GA_Print(ga);
  VecAssemblyBegin(vec);
  VecAssemblyEnd(vec);
  
  int ld;
  double *ptr;
  NGA_Access(ga,lo,hi,&ptr,&ld);
  PetscReal **d;
  int c=0;
  ierr = DAVecGetArray(da,vec,&d); CHKERRQ(ierr);
  for (int j = info.ys; j < info.ys+info.ym; ++j)
  {
    for (int i = info.xs; i < info.xs+info.xm; ++i)
    {
      if( d[j][i] != ptr[(i-info.xs)+ld*(j-info.ys)] )
        GA_Error("DA array is not equal to GA array",1);
//      printf("%d (%d,%d):\t%3.0f\t%3.0f\n", c, i, j, d[j][i], ptr[(i-info.xs)+ld*(j-info.ys)]);
      c++;
    }
  }
  ierr = DAVecRestoreArray(da,vec,&d); CHKERRQ(ierr);
  
  c=0;
  PetscReal *v;
  int start, end;
  VecGetOwnershipRange(vec, &start, &end);
  VecGetArray( vec, &v );
  for( int i = start; i < end; i++)
  {
//    printf("%d:\t%3.0f\t%3.0f\t%s\n", start, v[i-start], ptr[i-start], (v[i-start]-ptr[i-start]==0?"":"NO") );
  }
  VecRestoreArray( vec, &v );
  
  NGA_Release_update(ga,lo,hi);

  Vec gada;
//  VecCreateMPIWithArray(((PetscObject)da)->comm,da->Nlocal,PETSC_DETERMINE,ptr,&gada);
  VecView(gada,PETSC_VIEWER_STDOUT_SELF);
  
  GA_Destroy(ga);
  
  ierr = VecDestroy(vec); CHKERRQ(ierr);
  ierr = DADestroy(da); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

int main(int argc, char **args)
{
  PetscErrorCode  ierr;
  ierr = PetscInitialize(&argc, &args, (char *) 0, ""); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "Start %s\n", __FILE__); CHKERRQ(ierr);
  GA_Initialize();
  int stack = 1000;
  int heap = 1000;
  if( !MA_init(C_DBL, stack, heap))
      GA_Error("MA_init failed",stack+heap);
  
  size_t len = 1024;
  char name[1024];
  PetscGetHostName(name,len);
  ierr = PetscSynchronizedPrintf(PETSC_COMM_WORLD,"Hello from %s\n", name); CHKERRQ(ierr);
  ierr = PetscSynchronizedFlush(PETSC_COMM_WORLD); CHKERRQ(ierr);
  
  ierr = vizGA2DA(); CHKERRQ(ierr);
  
  GA_Terminate();
  ierr = PetscPrintf(PETSC_COMM_WORLD, "End %s\n", __FILE__); CHKERRQ(ierr);
  ierr = PetscFinalize(); CHKERRQ(ierr);
}
