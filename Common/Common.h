#ifndef COMMON_H_
#define COMMON_H_

// ignore warnings in other project headers
#pragma GCC diagnostic ignored "-Wpragmas"
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wcast-qual"
#include "ga.h"
#include "macdecls.h"
#include "petsc.h"
#pragma GCC diagnostic pop

#define LINE() { \
  int rank; \
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank); \
  printf("\n%s[%d]%s: %d\n", rank==1?"\t":"",rank,__FILE__, __LINE__); \
}

// Arbitrary tuple of at most three real numbers
// in 2D the 'z' coor is ignored
typedef struct {
  PetscReal x, y, z;
} Coor;

typedef struct {
  PetscInt x, y, z;
} iCoor;

typedef struct _AABB {
  Coor lo;
  Coor hi;
} AABB;

// Coordinate <--> index
inline void CoorToIndex( const Coor origin, const Coor dh, const Coor point, iCoor *P );
inline void CoorToIndex2(const Coor origin, const Coor dh, const Coor point, iCoor *P, Coor* p );
inline void IndexToCoor( const Coor origin, const Coor dh, const iCoor index, Coor *p );
inline PetscBool AABBPointInBox( const AABB box, const Coor p );

static const PetscReal Tensor1[][3] = {
  {0.5, 0.0, 0.0}, // u
  {0.0, 0.5, 0.0}, // v
  {0.0, 0.0, 0.5}  // w
};

static const PetscReal Tensor2[][3] = {
  {0.0,0.0,0.0}, // xx
  {0.5,0.5,0.0}, // xy
  {0.5,0.0,0.5}, // xz
  {0.0,0.0,0.0}, // yy
  {0.0,0.5,0.5}, // yz
  {0.0,0.0,0.0}  // zz
};

typedef enum {CELL_CENTER, U_FACE, V_FACE, W_FACE} VelFace;
static const iCoor STAGGERED_GRID[4] = {
    { 0, 0, 0},
    {-1, 0, 0},
    { 0,-1, 0},
    { 0, 0,-1}
};

static const int STAR[][3] = { {1,0,0},{-1, 0, 0},
                               {0,1,0},{ 0,-1, 0},
                               {0,0,1},{ 0, 0,-1} };

#define DCellInit() DCellInitialize(&argc,&args, __FILE__)
PetscErrorCode DCellInitialize(int *argc,char ***args, const char sourcefile[]);
PetscErrorCode DCellFinalize(void);
PetscErrorCode DCellGetWorkingDirectory( char* workingdirectory );

// Start offset of first file, so temporal files are alphabetically sorted properly
static const int FILE_COUNT_START = 10000;

PetscErrorCode VecWrite( Vec vec, const char *name, int t );
PetscErrorCode MatWrite( Mat mat, const char *name, int t );

typedef struct _LeastSq *LeastSq;
PetscErrorCode LeastSqCreate( int Np, PetscBool is2D, LeastSq *ls );
PetscErrorCode LeastSqDestroy( LeastSq ls );
PetscErrorCode LeastSqGetVecs( LeastSq ls, double **s, double **r, double **g, int *len );
PetscErrorCode LeastSqSolve( LeastSq ls );
PetscErrorCode LeastSqSetNumPoints( LeastSq ls, int n );

PetscErrorCode InterpolateVelocity2D( const int udof, PetscReal  ***field, const Coor X, Coor *vel );
PetscErrorCode InterpolateVelocity3D( const int udof, PetscReal ****field, const Coor X, Coor *vel );

// Array
typedef struct _Array *Array;
PetscErrorCode ArrayCreate( const char name[], int elemSize, Array *array );
PetscErrorCode ArrayDestroy( Array a );
PetscErrorCode ArraySetName( Array a, const char name[] );
PetscErrorCode ArraySetSize( Array a, int size );
PetscErrorCode ArrayAppend( Array a, void *elem );
PetscErrorCode ArrayAppendPtr( Array a, void *elem );
PetscErrorCode ArrayGet( Array a, int i, void *elem );
PetscErrorCode ArrayGetP( Array a, int i, void *elem );
PetscErrorCode ArrayWrite( Array a, int t );
PetscErrorCode ArrayZero( Array a );
PetscErrorCode ArrayDelete1( Array a, int idx );
PetscErrorCode ArrayCopy( Array src, Array copy );
PetscErrorCode ArrayDuplicate( Array a, Array *newarray );
PetscErrorCode ArrayGetCoor( Array a, iCoor pos, void *elem);
PetscErrorCode ArrayGetCoorP( Array a, iCoor pos, void *elem);
PetscErrorCode ArraySetCoor( Array a, iCoor shift, iCoor size );
int   ArrayLength( Array a );
void* ArrayGetData( Array a );

// Global Arrays - PETSc integration
typedef int GA; // todo: use GA type instead of int
PetscErrorCode GACreate( DM da, int *ga );
PetscErrorCode GAGetVec(int ga, Vec vec ); // GA --> PETSc
PetscErrorCode GAPutVec(Vec vec, int ga ); // PETSc --> GA
PetscErrorCode GAGather(int ga, Array c, Array v); // GA --> val[]
PetscErrorCode GAScatterAcc(int ga, Array c, Array v);// val[] +-> GA

// MemCache
typedef struct _MemCache *MemCache;
PetscErrorCode MemCacheCreate( const char name[], size_t elemsize, size_t chunksize, MemCache *mc );
PetscErrorCode MemCacheDestroy( MemCache mc );
PetscErrorCode MemCacheAlloc( MemCache mc, void *elem );
PetscErrorCode MemCacheFree( MemCache mc, void *elem );

// Heap
typedef struct _Heap *Heap;
typedef PetscReal (*Comparator)(void* parent, void* child);
typedef void (*HeapPrintNode)(void* node);
PetscErrorCode HeapCreate( const char name[], Comparator cmp, Heap *heap);
PetscErrorCode HeapDestroy( Heap h );
PetscErrorCode HeapClear( Heap h, MemCache mc );
PetscErrorCode HeapInsert( Heap h, void *item );
PetscErrorCode HeapPop( Heap h, void *topItem );
PetscErrorCode HeapPeek( Heap h, void *topItem );
PetscErrorCode HeapCheck( Heap h );
PetscErrorCode HeapPrint( Heap h, HeapPrintNode printNode );
int HeapSize( Heap h );
PetscBool HeapIsEmpty( Heap h );
PetscErrorCode  HeapRegisterEvents(void);

// UniqueID
typedef struct _UniqueID *UniqueID;
typedef int UniqueIDType;
PetscErrorCode UniqueIDCreate( UniqueID *uid );
PetscErrorCode UniqueIDDestroy( UniqueID uid );
PetscErrorCode UniqueIDSetStartCount( UniqueID uid, UniqueIDType maxID );
PetscErrorCode UniqueIDGenerate( UniqueID uid, UniqueIDType *id );

// SpatialIndex
typedef struct _SpatialIndex *SpatialIndex;
PetscErrorCode SpatialIndexCreate( const char name[], SpatialIndex *sidx );
PetscErrorCode SpatialIndexSetDomain( SpatialIndex sidx, Coor lo, Coor hi, Coor dh );
PetscErrorCode SpatialIndexDestroy( SpatialIndex sidx );
PetscErrorCode SpatialIndexInsertPoint( SpatialIndex sidx, Coor pt, void *data );
PetscErrorCode SpatialIndexInsertBox( SpatialIndex sidx, AABB box, void *item );
PetscErrorCode SpatialIndexQueryPoints( SpatialIndex sidx, Coor center, PetscReal radius, const int MAXLEN, int *len, void *items[] );
PetscErrorCode SpatialIndexQueryPointsBox( SpatialIndex sidx, AABB box, Array *items );
PetscErrorCode SpatialIndexClear( SpatialIndex sidx );
PetscErrorCode SpatialIndexInsertBox( SpatialIndex sidx, AABB box, void *item );
PetscErrorCode SpatialIndexCollide( SpatialIndex sidx, AABB box, void *items );
PetscErrorCode SpatialIndexPrint( SpatialIndex sidx );

#endif /* COMMON_H_ */
