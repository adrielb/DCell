#include "MicrofluidicSimulator.h"

void Divergence( UserContext *uc, Vec u, Vec v, Vec Div);
void Gradient( UserContext *uc, Vec P, Vec Px, Vec Py );
PetscErrorCode AssembleMatrices( UserContext *uc, Mat A, Mat B, Mat P );

#undef __FUNCT__
#define __FUNCT__ "PressureIncrement"
PetscInt EVENT_PressureIncrement;
PetscErrorCode PressureIncrement(UserContext *uc, Vec u, Vec v)
{
  PetscInt numVecsAllocated=10;
  Vec us, vs, div, p, px, py, phi, phiX, phiY, rhs, *vecs;
  Mat matVel, matVelRHS, matPhi;
  KSP kspVel, kspPhi;
  PetscReal dt;
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  PetscLogEventBegin(EVENT_PressureIncrement,0,0,0,0);
  
  VecDuplicateVecs(u,numVecsAllocated,&vecs);
  us=vecs[0]; vs=vecs[1]; div=vecs[2]; p=vecs[3]; px=vecs[4];  py=vecs[5];
  phi=vecs[6]; phiX=vecs[7]; phiY=vecs[8]; rhs=vecs[9];

  uc->DELTA_T = uc->DELTA_X * uc->DELTA_X * uc->DENSITY / (4 * uc->VISCOSITY);
  dt = uc->DELTA_T;
    
  MatCreate(PETSC_COMM_SELF, &matVel);
  MatCreate(PETSC_COMM_SELF, &matVelRHS);
  MatCreate(PETSC_COMM_SELF, &matPhi);
  MatSetSizes(matVel,    uc->numNodes, uc->numNodes, PETSC_DECIDE, PETSC_DECIDE);
  MatSetSizes(matVelRHS, uc->numNodes, uc->numNodes, PETSC_DECIDE, PETSC_DECIDE);
  MatSetSizes(matPhi,    uc->numNodes, uc->numNodes, PETSC_DECIDE, PETSC_DECIDE);
  MatSetType(matVel,    MATUMFPACK);
  MatSetType(matVelRHS, MATSEQAIJ);
//  MatSetType(matPhi,    MATUMFPACK);
  MatSetType(matPhi, MATSEQAIJ);
  MatSeqAIJSetPreallocation(matVel,   5,PETSC_NULL);
  MatSeqAIJSetPreallocation(matVelRHS,5,PETSC_NULL);
  MatSeqAIJSetPreallocation(matPhi,   5,PETSC_NULL);
  AssembleMatrices( uc, matVel, matVelRHS, matPhi);
  
  int s = 256;
  char tmp_dir[256], filename[256];
  PetscViewer view;
  ierr = PetscGetTmp(PETSC_COMM_SELF,tmp_dir,s); CHKERRQ(ierr);
  printf("tmp:%s\n",tmp_dir);
  sprintf(filename,"%s/matVel.dat",tmp_dir);
  ierr = PetscViewerASCIIOpen(PETSC_COMM_SELF, filename, &view); CHKERRQ(ierr);
  ierr = MatView(matPhi, view); CHKERRQ(ierr);
  ierr = PetscViewerDestroy(view); CHKERRQ(ierr);
  
  PC pc;
  KSPCreate(PETSC_COMM_SELF, &kspVel);
  KSPSetType(kspVel, KSPPREONLY);
  KSPGetPC(kspVel, &pc);
  PCSetType(pc, PCLU);
  KSPSetOperators(kspVel, matVel, matVel, SAME_PRECONDITIONER);
//  KSPSetInitialGuessNonzero(kspVel, PETSC_TRUE);
  KSPSetUp(kspVel);
  KSPCreate(PETSC_COMM_SELF, &kspPhi);
  KSPSetType(kspPhi, KSPPREONLY);
  KSPGetPC(kspPhi, &pc);
  PCSetType(pc, PCLU);
  KSPSetOperators(kspPhi, matPhi, matPhi, SAME_PRECONDITIONER);
//  KSPSetInitialGuessNonzero(kspPhi, PETSC_TRUE);
  KSPSetUp(kspPhi);

  for( int i = 0; i < uc->numBC; i++ )
  {
    VecSetValues(p,1,&uc->bcNodes[i].nodeIndex,&uc->bcNodes[i].pressureBC,INSERT_VALUES);
  }
  
  for( int i = 0; i < 20; i++)
  {
//  grad p = [ px, py ]
    Gradient(uc, p, px, py);
WriteResult(uc,p ,"p");
WriteResult(uc,px,"px");
WriteResult(uc,py,"py");

//  Solve [q - Laplace] . us = [q + Laplace] . un
    MatMult(matVelRHS, u, rhs);
    KSPSolve( kspVel, rhs, us);
//  Solve [q - Laplace] . vs = [q + Laplace] . vn
    MatMult(matVelRHS, v, rhs);
    KSPSolve( kspVel, rhs, vs);
WriteResult(uc,us,"us");
WriteResult(uc,vs,"vs");
    
//  div us = us:x + vs:y
    Divergence( uc, us, vs, div);
double norm;
VecNorm(div,NORM_2,&norm);
printf("[%d] %f\n", i, norm);
//  Solve Laplace(phi) = rho/dt * div 
    VecScale( div, uc->DENSITY / dt );
//  Pressure BC
    for( int i = 0; i < uc->numBC; i++ )
    {
      VecSetValues(div,1,&uc->bcNodes[i].nodeIndex,&uc->bcNodes[i].pressureBC,INSERT_VALUES);
    }
WriteResult(uc,div,"div");
    KSPSolve( kspPhi, div, phi);
WriteResult(uc,phi,"phi");

//  grad phi = [ phiX, phiY ]
    Gradient(uc, phi, phiX, phiY);
WriteResult(uc,phiX,"phiX");
WriteResult(uc,phiY,"phiY");
//  Velocity Correction
//  u = us - dt/rho * phiX
    VecWAXPY( u, -dt/uc->DENSITY, phiX, us);
    VecWAXPY( v, -dt/uc->DENSITY, phiY, vs);
WriteResult(uc,u,"u");
WriteResult(uc,v,"v");
//  Pressure Correction
//  div = dt/rho * div  to cancel previous scaling
//  p = -mu / 2 * div + phi
    VecWAXPY( p, -dt/uc->DENSITY * uc->VISCOSITY / 2., div, phi);
//WriteResult(uc,p,"p");
Divergence( uc, u, v, div);
VecNorm(div,NORM_2,&norm);
printf("[%d] %f\n", i, norm);
  }
  
  VecDestroyVecs(vecs, numVecsAllocated);
  MatDestroy(matVel);
  MatDestroy(matVelRHS);
  MatDestroy(matPhi);
  KSPDestroy(kspPhi);
  KSPDestroy(kspVel);
  PetscLogEventEnd(EVENT_PressureIncrement,0,0,0,0);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "AssembleMatrices"
PetscInt EVENT_AssembleMatrices;
PetscErrorCode AssembleMatrices( UserContext *uc, Mat A, Mat B, Mat P )
{
  int i, j;
  Node *n;
  PetscReal q = (2 * uc->DENSITY) / ( uc->VISCOSITY * uc->DELTA_T );
  PetscReal dx = uc->DELTA_X, h = 1/(dx*dx);
  PetscReal valH[5] = {-h,-h,-h,-h,q+4*h};
  PetscReal valR[5] = { h, h, h, h,q-4*h};
  PetscReal valN[5] = { h, h, h, h, -4*h};
  PetscReal valD[5] = { 0, 0, 0, 0,    1};
  PetscErrorCode ierr;
// 2 / (mu dt) u* - laplace(u*)  
  PetscFunctionBegin;
  PetscLogEventBegin(EVENT_AssembleMatrices,0,0,0,0);
  
  printf("DENSITY\t%f\n",uc->DENSITY);
  printf("VISCO\t%f\n",uc->VISCOSITY);
  printf("DELTA_T\t%f\n",uc->DELTA_T);
  printf("DELTA_X\t%f\n",uc->DELTA_X);
  printf("h\t\t%f\n",h);
  printf("q\t\t%f\n",q);
  
  for( i = 0; i < uc->numNodes; i++ )
  {
    n = &uc->nodes[i];
    if( n->isInterior == PETSC_TRUE ) { 
//    finite differences for interior nodes
      ierr = MatSetValues(A, 1, &i, 5, n->star, valH, INSERT_VALUES); CHKERRQ(ierr);
      ierr = MatSetValues(B, 1, &i, 5, n->star, valR, INSERT_VALUES); CHKERRQ(ierr);
      ierr = MatSetValues(P, 1, &i, 5, n->star, valN, INSERT_VALUES); CHKERRQ(ierr);
    }
  }
  
  for( i = 0; i < uc->numNodes; i++ )
  {
    n = &uc->nodes[i];
    if( n->isInterior != PETSC_TRUE ) {
//    no-slip condition for boundary velocity nodes ( v = u = 0 )
      ierr = MatSetValues(A, 5, n->star, 1, &i, valD, INSERT_VALUES); CHKERRQ(ierr);
      ierr = MatSetValues(B, 5, n->star, 1, &i, valD, INSERT_VALUES); CHKERRQ(ierr);
//    n . del( phi ) = 0
      valN[4] = -1*n->numNei*h;
      ierr = MatSetValues(P, 1, &i, 5, n->star, valN, INSERT_VALUES); CHKERRQ(ierr);
    }
  }
  
  int idx;
  for (i = 0; i < uc->numBC; ++i)
  {
    idx = uc->bcNodes[i].nodeIndex;
    n = &uc->nodes[idx];
//  phi = p at pressure BC
    ierr = MatSetValues(P,1,&idx,5,n->star,valD,INSERT_VALUES); CHKERRQ(ierr);
//    ierr = MatSetValues(P,5,n->star,1,&idx,valD,INSERT_VALUES); CHKERRQ(ierr);
  }   
  
  MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
  MatAssemblyBegin(B, MAT_FINAL_ASSEMBLY);
  MatAssemblyBegin(P, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(B, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(P, MAT_FINAL_ASSEMBLY);
  PetscLogEventEnd(EVENT_AssembleMatrices,0,0,0,0);
  PetscFunctionReturn(0);
}

void Gradient( UserContext *uc, Vec P, Vec Px, Vec Py )
{
  PetscReal *p, *px, *py;
  PetscInt N, S, E, W, C;
  Node *n;
  
  VecGetArray(P, &p);
  VecGetArray(Px, &px);
  VecGetArray(Py, &py);
  
  for( int i = 0; i < uc->numNodes; i++ )
  {
    n = &uc->nodes[i];
    
    S = n->star[0];
    W = n->star[1];
    E = n->star[2];
    N = n->star[3];
    C = n->star[4];
    
    if( n->isInterior) {
      px[i] = (p[E] - p[W]) / (2 * uc->DELTA_X);
      py[i] = (p[N] - p[S]) / (2 * uc->DELTA_X);
    } else {
      px[i] = 0;
      py[i] = 0;
    }
  }
  
  VecRestoreArray(P,&p);
  VecRestoreArray(Px, &px);
  VecRestoreArray(Py, &py);
}

void Divergence( UserContext *uc, Vec U, Vec V, Vec Div) // ux + vy
{
  int i;
  Node *n;
  PetscInt E,W,N,S;
  PetscReal uE, uW, vN, vS, *u, *v, *div;
  
  VecGetArray( U, &u);
  VecGetArray( V, &v);
  VecGetArray( Div, &div);
  for( i = 0; i < uc->numNodes; i++ )
  {
    n = &uc->nodes[i];
    S = n->star[0];
    W = n->star[1];
    E = n->star[2];
    N = n->star[3];
    
    if( n->isInterior )
      div[i] = ( u[E] - u[W] + v[N] - v[S] ) / (2 * uc->DELTA_X);
    else 
      div[i] = 0;
  }
  
  VecRestoreArray( Div, &div);
  VecRestoreArray( V, &v);
  VecRestoreArray( U, &u);
}

void Register_PressureIncrement()
{
  PetscLogEventRegister(&EVENT_PressureIncrement,"PressureIncrement", 0);
  PetscLogEventRegister(&EVENT_AssembleMatrices,"AssembleMatrices", 0);
}