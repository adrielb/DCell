#ifndef DCELL_H_
#define DCELL_H_
#include "LevelSetMethod.h"
#include "Reaction.h"
#include "petscda.h"

//TODO: try to make 2D and 3D coalesce
//TODO: if global coordinates are int R3 (instead of Z3), then 
typedef struct _DCell {
  LevelSet lsPlasmaMembrane;
//  LevelSet2D lsNuclearMembrane;
  Grid u, v, w; // Advection solution vectors
  Reaction rxn; // Reaction fuctions/jacobian
} *DCell;

PetscErrorCode DCellCreate(Reaction rxn, DCell *cell);

/// ??? /// TODO what is this?
void DCellAssembleDiffusion( LevelSet ls, Mat mat );
void AssembleDiffusion2DBC( Mat mat, MatStencil row, int* lens, 
    PetscReal D, PetscReal dx2, PetscReal dy2 );

#endif /*DCELL_H_*/
