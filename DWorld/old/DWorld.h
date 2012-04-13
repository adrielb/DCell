#ifndef DWORLD_H_
#define DWORLD_H_

#include "Reaction.h"
#include "Utilities.h"
#include "petscda.h"
#include "petscts.h"
#include "FluidField.h"

typedef struct _MassField {
  DA da; // this DA has multiple DOF for each chemical specie
  Coor len;  // Length of each dimension (typically microns)
  Coor d;    // Grid spacing [dx, dy, dz]
  Vec conc;  // Concentration  
} *MassField;

typedef struct _DWorld { //Application context
  FluidField fluid;
  MassField mass;
  Mat L, J; // the linear and non-linear part
  Array cells; // cells this processor owns
} *DWorld;

PetscErrorCode DWorldCreate( Reaction rxn, iCoor num, Coor len, DWorld *world);
PetscErrorCode DWorldDestroy( DWorld world);

PetscErrorCode RHSFunction2D(TS ts,PetscReal t,Vec globalin,Vec globalout,void *ctx);
PetscErrorCode RHSJacobian2D(TS ts,PetscReal t,Vec x,Mat *AA,Mat *BB,MatStructure *flag,void *ptr);

PetscErrorCode RHSFunction(TS ts,PetscReal t,Vec globalin,Vec globalout,void *ctx);
PetscErrorCode RHSJacobian(TS ts,PetscReal t,Vec x,Mat *AA,Mat *BB,MatStructure *flag,void *ptr);

PetscErrorCode Monitor( TS ts,PetscInt step,PetscReal time,Vec global,void *ctx);

PetscErrorCode AssembleAdvectiveTransport( DALocalInfo info, Mat mat, Coor d, PetscReal **u, PetscReal **v );

#endif /*DWORLD_H_*/
