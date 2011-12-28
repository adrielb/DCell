#ifndef REACTION_H_
#define REACTION_H_
#include "petsc.h"

typedef void (*ReactionFunction)(PetscReal *c, PetscReal *F);

typedef struct _Reaction {
  PetscInt dof; // Number of chemical species
  PetscReal *D; // Diffusion coefficient (list length dof)
  PetscReal *F; // Output of reaction function evaluation (length dof)
  ReactionFunction ComputeFunction;     
  ReactionFunction ComputeJacobian;   
  PetscReal *jac; // Values for evaluated Jacobian (up to dof^2)
  PetscInt *rows; 
  PetscInt *cols;
  PetscInt jac_length; //Length of jac,rows,cols array
} *Reaction;

PetscErrorCode ReactionCreate( PetscInt num, Reaction *rxn );
PetscErrorCode ReactionDestroy( Reaction rxn );
void ReactionSetFunction( Reaction rxn, ReactionFunction func);
void ReactionSetJacobian( Reaction rxn, ReactionFunction jac);
void ReactionUpdateFunction( Reaction rxn, PetscReal *F );
void ReactionUpdateJacobian( Reaction rxn, PetscReal *jac );

void ReactionFunction_Null( PetscReal *c, PetscReal *F);
void ReactionCreate_Test( int dof, Reaction *rxn );
void ReactionCreate_Barkley( Reaction *rxn );

#endif /*REACTION_H_*/
