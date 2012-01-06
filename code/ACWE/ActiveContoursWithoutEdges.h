#ifndef ACTIVECONTOURSWITHOUTEDGES_H_
#define ACTIVECONTOURSWITHOUTEDGES_H_

#include "petscksp.h"
#include "LevelSetMethod.h"

typedef struct _ACWE {
  PetscReal dt;
  PetscReal mu, l1, l2, c1, c2;
  LevelSet  ls;
  Grid      vel;
  Array     idxArray;
} *ACWE;

PetscErrorCode ACWECreate( Grid img, ACWE *uctx );
PetscErrorCode ACWEDestroy( ACWE uc );

PetscErrorCode SingleStep( ACWE uc, Grid img );
PetscErrorCode IndexInterior(ACWE uc);
PetscErrorCode CalculateC1C2( Array idxArray, PetscReal *phi, PetscReal *image,
  PetscReal *c1, PetscReal *c2);



PetscErrorCode RegisterEvents_ACWE();

#endif /*ACTIVECONTOURSWITHOUTEDGES_H_*/
