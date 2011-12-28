#ifndef POISSONSOLVER_H_
#define POISSONSOLVER_H_

#include "petscksp.h"

PetscErrorCode Generate2DLaplacianPeriodicBC(  PetscInt d1, PetscInt d2, Mat *m );
PetscErrorCode Generate2DLapacian( PetscInt d1, PetscInt d2, Mat *mat );
PetscErrorCode GenerateLaplacian2DNoBC( PetscInt d1, PetscInt d2, Mat *m );
PetscErrorCode SpectralMethod( int nn, PetscLogDouble *time );
PetscErrorCode SpectralMethod2D();

#endif /*POISSONSOLVER_H_*/
