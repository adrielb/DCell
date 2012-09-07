#ifndef IMMERSEDINTERFACEMETHOD_H_
#define IMMERSEDINTERFACEMETHOD_H_

#include "LevelSetMethod.h"

typedef struct _IIMIrregularNode IIMIrregularNode;
typedef void (*InterfacialForce)( IIMIrregularNode *n, void *context);
void InterfacialForceSurfaceTension( IIMIrregularNode *n, void *context );
void InterfacialForceElastic( IIMIrregularNode *n, void *context );
void InterfacialForceSimple( IIMIrregularNode *n, void *context );


typedef struct _IIM *IIM;
PetscErrorCode IIMCreate( PetscBool is2D, IIM *iim );
PetscErrorCode IIMDestroy( IIM iim );
PetscErrorCode IIMUpdateRHS( IIM iim, LevelSet ls, int ga );
PetscErrorCode IIMCorrectVelocity( IIM iim, const Coor X, Coor *vel );
PetscErrorCode IIMWriteIrregularNodeList( IIM iim, char *name, int idx );

PetscErrorCode IIMSetForceComponents( IIM iim, InterfacialForce F );
PetscErrorCode IIMSetForceContext( IIM iim, void *context );
PetscErrorCode IIMSetViscosity( IIM iim, PetscReal mu );
PetscErrorCode IIMSetEps( IIM iim, PetscReal eps );
PetscErrorCode IIMSetNp( IIM iim, int Np );
PetscErrorCode IIMSetFluidCoordinateSystem( IIM iim, Coor origin, Coor df );


typedef struct _LocalCoor *LocalCoor;
PetscErrorCode LocalCoorCreate( LocalCoor *lc );
PetscErrorCode LocalCoorDestroy( LocalCoor lc );
void LocalCoorSetLength( LocalCoor lc, PetscInt len);
void LocalCoorSolve( LocalCoor lc, Coor dh, IIMIrregularNode *N, IIMIrregularNode *nodes[] );
void LocalCoorGetVecs( LocalCoor lc, PetscReal **s, PetscReal **n, PetscReal **r );
void LocalCoor3DTangential( IIMIrregularNode *n );

#endif /*IMMERSEDINTERFACEMETHOD_H_*/
