#ifndef STOKESFLOW_H_
#define STOKESFLOW_H_

#include "petscksp.h"
#include "Utilities.h"

struct node_type {
  PetscInt numNei;
  PetscInt nei[5];
  PetscInt star[5];
  PetscInt imageIndex;
  PetscTruth isInterior;
};

typedef struct node_type Node;

struct bcnode_type {
  PetscScalar pressureBC;
  PetscScalar concentrationBC;
  PetscInt    nodeIndex;
};

typedef struct bcnode_type BCNode;

typedef struct {
  char            *imageFileName;
  PetscInt        n, numrows, numcols;
  unsigned char   *filedata;
  BCNode          BCLabels[256];
  PetscInt        *imageToNode;
  PetscInt        numBC;
  BCNode          *bcNodes;
//  PetscInt        *bcToImage;
  PetscInt        numNodes;
  Node            *nodes;
  PetscReal       *imageResult;
  PetscInt        BACKGROUND_COLOR, 
                  FLUIDIC_LAYER_COLOR, 
                  SCALE_COLOR;
  PetscReal       VISCOSITY, DIFFUSION, 
                  DELTA_X,   DELTA_Z,
                  DENSITY,   GRAVITY,
                  SCALE_BAR, DELTA_T;
} UserContext;

/* Microfluidic Simulator
 */
PetscErrorCode ReadFile(char *name, PetscInt len, unsigned char **filedata);
PetscErrorCode InterpretOptions(UserContext *uc);
PetscErrorCode IndexFreeNodes( UserContext* ); 
PetscErrorCode DetermineScale( UserContext *uc );
PetscErrorCode Histogram( UserContext *uc );
PetscErrorCode WriteResult( UserContext *uc, Vec v, char *filename );
PetscErrorCode DestroyContext( UserContext* );
void MicrofluidicSimulator_RegisterEvents();

PetscErrorCode PressureIncrement(UserContext *uc, Vec u, Vec v);

/*  Diffusion Advection
 */
PetscErrorCode AssembleSSConcentration(UserContext *uc);
PetscErrorCode SolveConcentration( UserContext *uc );
PetscErrorCode AssembleDiffusionMatrix( UserContext *uc );
PetscErrorCode AssembleConcentrationRHS( UserContext *uc );

#endif /*STOKESFLOW_H_*/
