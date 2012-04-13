#include "Reaction.h"

#undef __FUNCT__
#define __FUNCT__ "ReactionCreate"
PetscInt EVENT_ReactionCreate;
PetscErrorCode ReactionCreate( int dof, Reaction *rxn)
{
  struct _Reaction *r;
  int memsize = dof*sizeof(PetscReal);
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  PetscLogEventBegin(EVENT_ReactionCreate,0,0,0,0);
//  PetscLogEventRegister(&EVENT_ReactionCreate,"ReactionCreate", 0);
  
  ierr = PetscNew(struct _Reaction, &r); CHKERRQ(ierr);
  ierr = PetscMalloc(memsize, &r->D); CHKERRQ(ierr);
  ierr = PetscMalloc(memsize, &r->F); CHKERRQ(ierr);
  ierr = PetscMalloc(dof*memsize, &r->jac); CHKERRQ(ierr);
  ierr = PetscMalloc(dof*dof*sizeof(int), &r->rows); CHKERRQ(ierr);
  ierr = PetscMalloc(dof*dof*sizeof(int), &r->cols); CHKERRQ(ierr);
  
  // Set default reactions to functions without kinetics
  r->ComputeFunction = ReactionFunction_Null;
  r->ComputeJacobian = ReactionFunction_Null;
  
  PetscMemzero(r->D,memsize);
  PetscMemzero(r->F,memsize);
  PetscMemzero(r->jac,memsize);

  r->dof = dof;
  r->jac_length = 0;
  
  *rxn = r; // Return contructed object
  
  PetscLogEventEnd(EVENT_ReactionCreate,0,0,0,0);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ReactionDestroy"
PetscErrorCode ReactionDestroy( Reaction rxn )
{
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  ierr = PetscFree(rxn->D); CHKERRQ(ierr);
  ierr = PetscFree(rxn->F); CHKERRQ(ierr);
  ierr = PetscFree(rxn->jac); CHKERRQ(ierr);
  ierr = PetscFree(rxn->rows); CHKERRQ(ierr);
  ierr = PetscFree(rxn->cols); CHKERRQ(ierr);
  ierr = PetscFree(rxn); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

void ReactionSetFunction( Reaction rxn, ReactionFunction func)
{
  rxn->ComputeFunction = func;
}

void ReactionSetJacobian( Reaction rxn, ReactionFunction jac)
{
  rxn->ComputeJacobian = jac;
}

void ReactionUpdateFunction( Reaction rxn, PetscReal *chem )
{
  rxn->ComputeFunction( chem, rxn->F );
}

void ReactionUpdateJacobian( Reaction rxn, PetscReal *chem )
{
  rxn->ComputeJacobian( chem, rxn->jac );
}

void ReactionFunction_Null( PetscReal *c, PetscReal *F)
{
  // Do nothing
}

//https://computation.llnl.gov/iscr/annual_report/fy2006/Content/subcontracts/holst.html

void ReactionCreate_Test( int dof, Reaction *rxn )
{
  Reaction r;
  ReactionCreate( dof, &r );
  
  for (int i = 0; i < dof; ++i) {
    r->D[i] = i+1;
  }
  
  r->jac_length = dof*dof;
  
  int count = 0;
  for (int i = 0; i < dof; ++i) {
    for (int j = 0; j < dof; ++j) {
      r->rows[count] = i;
      r->cols[count] = j;
      r->jac[count] = 0; //count+1;
      count++;
    }
  }
  
  *rxn = r;
}
void ReactionFunction_SimpleDegradation( PetscReal *c, PetscReal *F );
void ReactionCreate_SimpleDegradation( Reaction *rxn )
{
  int dof = 5;
  Reaction r;
  ReactionCreate( dof, &r );
  ReactionSetFunction(r, ReactionFunction_SimpleDegradation);
  
  for (int i = 0; i < dof; ++i) {
    r->D[i] = i;
  }
  
  r->jac_length = dof;
  
  for (int j = 0; j < dof; ++j) {
    r->rows[j] = j;
    r->cols[j] = j;
    r->jac[j] = -(j+1);
  }
  
  *rxn = r;
}

void ReactionFunction_SimpleDegradation( PetscReal *c, PetscReal *F )
{
  int dof = 5;
  for (int j = 0; j < dof; ++j) {
    F[j] = -j * c[j];
  }
}

void ReactionFunction_GlycolyticOscillator( PetscReal *c, PetscReal *F );
void ReactionJacobian_GlycolyticOscillator( PetscReal *c, PetscReal *F );
void ReactionCreate_GlycolyticOscillator( Reaction *rxn )
{
  Reaction r;
  
  ReactionCreate(2, &r);
  ReactionSetFunction(r, ReactionFunction_GlycolyticOscillator);
  ReactionSetJacobian(r, ReactionJacobian_GlycolyticOscillator);

  int count = 0, i, j;
  for (i = 0; i < 2; ++i) {
    for (j = 0; j < 2; ++j) {
      r->rows[count] = i;
      r->cols[count] = j;
      count++;
    }
  }
  
  *rxn = r;
}

void ReactionFunction_GlycolyticOscillator( PetscReal *c, PetscReal *F )
{
  const PetscReal a = 0.04, b = 0.4,
    x = c[0], y = c[1]; 
  
  F[0] = -x + a*y + x*x*y;
  F[1] =  b - a*y - x*x*y;
}

void ReactionJacobian_GlycolyticOscillator( PetscReal *c, PetscReal *j )
{
  const PetscReal a = 0.04, b = 0.4,
      x = c[0], y = c[1]; 
  
  j[0] = -1 + 2*x*y;
  j[1] = a + x*x;
  j[2] = -2*x*y;
  j[3] = -a - x*x;
}

void ReactionFunction_GrayScott( PetscReal *c, PetscReal *F);
void ReactionJacobian_GrayScott( PetscReal *c, PetscReal *j);
void ReactionCreate_GrayScott( Reaction *rxn)
{
  Reaction r;
  
  ReactionCreate(2,&r);
  ReactionSetFunction(r,ReactionFunction_GrayScott);
  ReactionSetJacobian(r,ReactionJacobian_GrayScott);
  
  r->D[0] = 8e-5;
  r->D[1] = 4e-5;
  
  int count = 0, i, j;
  for (i = 0; i < 2; ++i) {
    for (j = 0; j < 2; ++j) {
      r->rows[count] = i;
      r->cols[count] = j;
      count++;
    }
  }
  
  *rxn = r;
}

void ReactionFunction_GrayScott( PetscReal *c, PetscReal *F)
{
  const PetscReal g = 0.024, k = 0.06,
    u = c[0], v = c[1];
  
  F[0] = -u*v*v + g*(1-u);
  F[1] =  u*v*v - (g+k)*v;
}
void ReactionJacobian_GrayScott( PetscReal *c, PetscReal *j)
{
  const PetscReal g = 0.024, k = 0.06,
      u = c[0], v = c[1];

  j[0] = -g - v*v;
  j[1] = -2*u*v;
  j[2] = v*v;
  j[3] = -g - k + 2 *u*v;
}

/* Barkley, 1991 */
void ReactionFunction_Barkley( PetscReal *c, PetscReal *F );
void ReactionJacobian_Barkley( PetscReal *c, PetscReal *val);

void ReactionCreate_Barkley(Reaction *rxn)
{
  Reaction r;
  ReactionCreate( 2, &r);
  ReactionSetFunction(r, ReactionFunction_Barkley);
  ReactionSetJacobian(r, ReactionJacobian_Barkley);
  r->D[0] = 1; 
  r->D[1] = 0;
  r->jac_length=4;
  PetscInt rows[4] = {0,0,1,1},
           cols[4] = {0,1,0,1};
  PetscMemcpy(r->rows, &rows, 4*sizeof(PetscInt));
  PetscMemcpy(r->cols, &cols, 4*sizeof(PetscInt));
  *rxn = r;
}

void ReactionFunction_Barkley( PetscReal *c, PetscReal *F )
{
  const PetscReal eps = 0.002, a = 0.25, b = 0.001,
        u = c[0], v = c[1];
  
  F[0] = 1 / eps * u * ( 1 - u ) * ( u - 1 / a * ( v + b)); 
  F[1] = u - v;
}

void ReactionJacobian_Barkley( PetscReal *c, PetscReal *val)
{
  const PetscReal eps = 0.002, a = 0.25, b = 0.001,
        u = c[0], v = c[1];
  
  val[0] = -((b - 2*b*u + a*u*(-2 + 3*u) + v - 2*u*v)/(a*eps)); 
  val[1] = ((-1 + u)*u)/(a*eps);
  val[2] = 1;
  val[3] = -1;
}

void ReactionFunction_Turing( PetscReal *c, PetscReal *F);
void ReactionJacobian_Turing( PetscReal *c, PetscReal *val);
void ReactionCreate_Turing( Reaction *rxn)
{
  Reaction r;
  ReactionCreate( 2, &r);
  ReactionSetFunction(r, ReactionFunction_Turing);
  ReactionSetJacobian(r, ReactionJacobian_Turing);
  r->D[0] = 0.25;
  r->D[1] = 0.0625;
  *rxn = r;
}

void ReactionFunction_Turing( PetscReal *c, PetscReal *F)
{
  const PetscReal beta = 12, a = c[0], b = c[1], s = 0.05;
  
  F[0] = s * (16 - a*b);
  F[1] = s * (a*b - b - beta);
}

void ReactionJacobian_Turing( PetscReal *c, PetscReal *val)
{
  const PetscReal beta = 12, a = c[0], b = c[1], s = 0.05;
  
  val[0] = -b * s;
  val[1] = -a * s;
  val[2] =  b * s;
  val[3] = (-1+a) * s;
}