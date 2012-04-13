#ifndef UTILITIES_H
#define UTILITIES_H

#include "petscts.h"
#include "petscda.h"

// Arbitrary tuple of at most three real numbers
// in 2D the 'z' coor is ignored
typedef struct {
  PetscReal x, y, z;
} Coor;
//TODO: rename type this to R3 and Z3, then 'typedef R3 Lengths'
typedef struct {
  PetscInt x, y, z;
} iCoor;


static const PetscReal Tensor2[][3] = {
  {0.0,0.0,0.0}, // xx
  {0.5,0.5,0.0}, // xy
  {0.5,0.0,0.5}, // xz
  {0.0,0.0,0.0}, // yy
  {0.0,0.5,0.5}, // yz
  {0.0,0.0,0.0}  // zz
};

static const PetscReal Tensor1[][3] = {
  {0.5,   0,   0}, // u
  {  0, 0.5,   0}, // v
  {  0,   0, 0.5}  // w
};

/* DomainSpecs
 * Different from DALocalInfo since DALocalInfo may include interior and boundary nodes
 * DomainSpecs just specifies only the interior indexes and lengths   
 *    ............
 *    ............
 * ---..********..
 * ly ..********.. 
 * ---..********..
 *    ............
 *    ............
 *      |--lx--|
 */
typedef struct _DomainSpecs {
  Coor len; // real length of domain (in units) 
  Coor d;   // grid spacing [dx, dy dz]
  iCoor s,e; // global start and end index for domain
  PetscTruth is2D;  // whether to loop over 'z' coordinate
} *DomainSpecs;

static const PetscReal negone=-1., zero=0., half = 0.5, one = 1., two = 2., four = 4.;

PetscErrorCode RegisterUtilityEvents();
PetscErrorCode WriteVector( char *name, Vec v );
PetscErrorCode WriteVectorN( char *name, int i, Vec v );
PetscErrorCode WriteVectorArray( char *name, PetscInt len, PetscReal *data );
PetscErrorCode MatWrite( char *name, Mat mat );

typedef struct _WriteVec* WriteVec;
PetscErrorCode WriteVecCreate( Vec v, char* name, WriteVec *write );
PetscErrorCode WriteVecDestroy ( WriteVec wv );
PetscErrorCode WriteVecToDisk  ( WriteVec wv );
           int WriteVecGetCount( WriteVec wv );

/*
#define negone -1.
#define zero 0.
#define one 1.
#define two 2.
#define three 3.
#define four 4.
*/

typedef struct _LeastSq *LeastSq;

PetscErrorCode LeastSqCreate( int Np, PetscTruth is2D, LeastSq *ls );
PetscErrorCode LeastSqDestroy( LeastSq ls );
PetscErrorCode LeastSqGetVecs( LeastSq ls, double **s, double **r, double **g, int *len );
PetscErrorCode LeastSqSolve( LeastSq ls );
PetscErrorCode LeastSqSetNumPoints( LeastSq ls, int n );

//TODO: find a better place for GlobalArray interface
PetscErrorCode DACreateGlobalArray( DA da, int *GA, Vec *g );
PetscErrorCode VecGetValuesGA( DA, Vec vec, int ga, PetscReal *val, int **idx, int len );
#endif /*  UTILITIES_H  */
