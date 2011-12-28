#include "FluidField.h"
#include "ga.h"
#include "macdecls.h"
#include PETSC_DIR/include/private/daimpl.h"

/*
 * TODO: change the communicator indexing instead of the GA malloc 
 * 
No problem, here is the code:

// the numbers of processors per direction are (int) x_procs, y_procs, z_procs respectively 
// (no parallelization in direction 'dir' means dir_procs = 1)

MPI_Comm NewComm;
int MPI_Rank, NewRank, x,y,z;

// get rank from MPI ordering:
MPI_Comm_rank(MPI_COMM_WORLD, &MPI_Rank);

// calculate coordinates of cpus in MPI ordering:
x = MPI_rank / (z_procs*y_procs);
y = (MPI_rank % (z_procs*y_procs)) / z_procs;
z = (MPI_rank % (z_procs*y_procs)) % z_procs;

// set new rank according to PETSc ordering:
NewRank = z*y_procs*x_procs + y*x_procs + x;

// create communicator with new ranks according to PETSc ordering:
MPI_Comm_split(PETSC_COMM_WORLD, 1, NewRank, &NewComm);

// override the default communicator (was MPI_COMM_WORLD as default)
PETSC_COMM_WORLD = NewComm;

I hope, this will be useful for some of you.

Ciao,
Rolf
*/

//TODO: use a better test between MPI and UNI processes
#ifdef MPIUNI_INTPTR // petsc is using MPI UNI, GA library not available

#undef __FUNCT__
#define __FUNCT__ "DACreateGlobalArray"
PetscErrorCode DACreateGlobalArray( DA da, int *GA, Vec *g )
{
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  
  ierr = DACreateGlobalVector(da, g); CHKERRQ(ierr);
  *GA = 0;
  
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "VecGetValuesGA"
PetscLogEvent EVENT_VecGetValuesGA;
PetscErrorCode VecGetValuesGA( DA da, Vec vec, int ga, PetscReal *val, int **idx, int len )
{
  DALocalInfo info;
  int i;
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  PetscLogEventBegin(EVENT_VecGetValuesGA,0,0,0,0);
//  PetscLogEventRegister("VecGetValues", 0, &EVENT_VecGetValues);
  ierr = DAGetLocalInfo(da,&info); CHKERRQ(ierr);
  if( info.dim == 2 )
  {
    PetscReal **array;
    ierr = DAVecGetArray(da, vec, &array); CHKERRQ(ierr);
    for( i = 0; i < len; ++i) {
      val[i] = array[idx[i][1]][idx[i][0]];
    }
  } else {
    PetscReal ***array;
    ierr = DAVecGetArray(da, vec, &array); CHKERRQ(ierr);
    for( i = 0; i < len; ++i) {
      val[i] = array[idx[i][2]][idx[i][1]][idx[i][0]];
    }    
  }
  
  PetscLogFlops( 0 );
  PetscLogEventEnd(EVENT_VecGetValuesGA,0,0,0,0);
  PetscFunctionReturn(0);
}

#else // petsc is using standard MPI, GA library available

#undef __FUNCT__
#define __FUNCT__ "VecGetValuesGA"
PetscLogEvent EVENT_VecGetValuesGA;
PetscErrorCode VecGetValuesGA( DA da, Vec vec, int ga, PetscReal *val, int **idx, int len )
{
  PetscErrorCode ierr;
  
  
  PetscFunctionBegin;
  PetscLogEventBegin(EVENT_VecGetValuesGA,0,0,0,0);
//  PetscLogEventRegister("VecGetValues", 0, &EVENT_VecGetValues);
  
  NGA_Gather(ga, val, idx, len);
  
  PetscLogFlops( 0 );
  PetscLogEventEnd(EVENT_VecGetValuesGA,0,0,0,0);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DACreateGlobalArray"
PetscErrorCode DACreateGlobalArray( DA da, int *GA, Vec *g )
{
  PetscInt ndim,M,N,P,m,n,p;
  const PetscInt *lx, *ly, *lz;
  int ga, *map;
  DALocalInfo info;
  int rank;
  int i;
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
  ierr = DAGetLocalInfo(da, &info); CHKERRQ(ierr);
  ierr = DAGetInfo(da,&ndim,&M,&N,&P,&m,&n,&p,0,0,0,0); CHKERRQ(ierr);
  ierr = DAGetOwnershipRanges(da, &lx, &ly, &lz); CHKERRQ(ierr);
  
  ga = GA_Create_handle();
  
  if( ndim == 2 )
  {
    int dims[2] = {N,M};
    int block[2] = {n,m};
    GA_Set_data(ga,2,dims,MT_DBL);
    PetscMalloc( sizeof(int)*(m+n), &map);
    map[0] = 0;
    for( i = 1; i < n; i++ )
    {
      map[i] = ly[i-1] + map[i-1];
    }
    map[n] = 0;
    for( i = n+1; i < m+n; i++ )
    {
      map[i] = lx[i-n-1] + map[i-1];
    }
    GA_Set_irreg_distr(ga,map,block);
    ierr = PetscFree(map); CHKERRQ(ierr);
    ierr = GA_Allocate( ga );
    
    int lo[2], hi[2];
    NGA_Distribution(ga,rank,lo,hi);
    if( lo[1] != info.xs || hi[1] != info.xs+info.xm-1 ||
        lo[0] != info.ys || hi[0] != info.ys+info.ym-1 )
    {
      PetscPrintf(PETSC_COMM_SELF,"[%d] lo:(%2d,%2d)  hi:(%2d,%2d) \t DA: (%2d,%2d), (%2d, %2d)\n",
          rank, lo[1], lo[0], hi[1], hi[0], info.xs, info.ys, info.xs+info.xm-1, info.ys+info.ym-1);
      GA_Error( "GA Distribution does not match DA distribution", 1);
    }
  }
  if( ndim == 3 )
  {
    int dims[3] = {P,N,M};
    int block[3] = {p,n,m};
    GA_Set_data(ga,ndim,dims,MT_DBL); // TODO: how to ensure MT_DBL is of same type as PetscReal?
    ierr = PetscMalloc( sizeof(int)*(m+n+p), &map); CHKERRQ(ierr);
    int i = 0;
    map[0] = 0;
    for( i = 1; i < p; i++ )  map[i] = map[i-1] + lz[i-1];
    map[p] = 0;
    for( i=p+1; i < n+p; i++ )  map[i] = map[i-1] + ly[i-1-p];
    map[n+p] = 0;
    for( i=n+p+1; i < m+n+p; i++ )  map[i] = map[i-1] + lx[i-1-p-n];
    GA_Set_irreg_distr(ga,map,block);
    ierr = PetscFree(map); CHKERRQ(ierr);
    ierr = GA_Allocate( ga );
    
    int lo[3], hi[3];
    NGA_Distribution(ga,rank,lo,hi);
    if( lo[2] != info.xs || hi[2] != info.xs+info.xm-1 ||
        lo[1] != info.ys || hi[1] != info.ys+info.ym-1 ||
        lo[0] != info.zs || hi[0] != info.zs+info.zm-1 )
    {
      PetscPrintf(PETSC_COMM_SELF,"[%d] lo:(%d,%d,%d)  hi:(%d,%d,%d) \t DA: (%d,%d,%d), (%d,%d,%d)\n",
          rank, lo[2],lo[1],lo[0], hi[2],hi[1],hi[0], 
          info.xs, info.ys, info.zs,
          info.xs+info.xm-1, info.ys+info.ym-1, info.zs+info.zm-1);
      GA_Error( "GA Distribution does not match DA distribution", 1);
    }
  }
  //TODO: better error handling for various GA failures
  if( !ierr ) GA_Error("\n\n\nga allocaltion failed\n\n",ierr);
  if( !ga ) GA_Error("\n\n\n ga null \n\n",1);
  if( rank != GA_Nodeid() ) GA_Error("MPI rank does not match GA_Nodeid()",1);
  
  PetscReal *ga_ptr;
  int lo[3], hi[3], ld;
  NGA_Distribution(ga,rank,lo,hi);
  NGA_Access(ga,lo,hi,&ga_ptr,&ld); //TODO: why does NGA_Access need lo,hi coordinates?
  /* Code below is taken directly from DACreateGlobalVector()
   * only change is VecCreateMPI() -> VecCreateMPIWithArray( ..., ga_ptr, g )
   */   
  PetscValidHeaderSpecific(da,DM_COOKIE,1);
  PetscValidPointer(g,2);
  ierr = VecCreateMPIWithArray(((PetscObject)da)->comm,da->Nlocal,PETSC_DETERMINE,ga_ptr,g); CHKERRQ(ierr);
  ierr = PetscObjectCompose((PetscObject)*g,"DA",(PetscObject)da);CHKERRQ(ierr);
  ierr = VecSetLocalToGlobalMapping(*g,da->ltogmap);CHKERRQ(ierr);
  ierr = VecSetLocalToGlobalMappingBlock(*g,da->ltogmapb);CHKERRQ(ierr);
  ierr = VecSetBlockSize(*g,da->w);CHKERRQ(ierr);
  ierr = VecSetOperation(*g,VECOP_VIEW,(void(*)(void))VecView_MPI_DA);CHKERRQ(ierr);
  ierr = VecSetOperation(*g,VECOP_LOADINTOVECTOR,(void(*)(void))VecLoadIntoVector_Binary_DA);CHKERRQ(ierr);

  /* Plan B if 'ga_ptr' invalid after call to NGA_Release
  ierr = DACreateGlobalVector(da,&vec); CHKERRQ(ierr);
  ierr = VecPlaceArray(vec,ga_ptr); CHKERRQ(ierr);
  */
//  NGA_Release(ga,lo,hi);
  
  *GA = ga;
  
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "VecGetGlobalArray"
PetscErrorCode VecGetGlobalArray( Vec vec, int ga )
{
  PetscReal *ga_ptr;
  int lo[3], hi[3], ld;
  int rank = GA_Nodeid();
  PetscErrorCode ierr;
  
  
  PetscFunctionBegin;
  
  NGA_Distribution(ga,rank,lo,hi);
  NGA_Access(ga,lo,hi,&ga_ptr,&ld); //TODO: why does NGA_Access need lo,hi coordinates?
  ierr = VecPlaceArray(vec,ga_ptr); CHKERRQ(ierr);
  
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "VecRestoreGlobalArray"
PetscErrorCode VecRestoreGlobalArray( Vec vec, int ga )
{
  int lo[3], hi[3];
  int rank = GA_Nodeid();
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  
  NGA_Distribution(ga,rank,lo,hi);
  NGA_Release_update(ga,lo,hi);
  ierr = VecResetArray(vec); CHKERRQ(ierr);
  
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "VecDestroyGlobalArray"
PetscErrorCode VecDestroyGlobalArray(Vec vec, int ga)
{
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  
  PetscFunctionReturn(0);
}

#endif // end MPIUNI test
