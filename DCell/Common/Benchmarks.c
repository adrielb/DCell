#include "Main.h"

PetscErrorCode RegisterEvents_Benchmarks(  );
PetscErrorCode MultipleHornerSchemes(int n, unsigned int iter);
PetscErrorCode RandomAccess( PetscInt s );
PetscErrorCode HornerScheme(unsigned int iter, double *X);

#undef __FUNCT__
#define __FUNCT__ "PetscMain"
PetscErrorCode PetscMain(  )
{
  PetscErrorCode ierr;
  int stage;
  double x;
  
  PetscFunctionBegin;
  RegisterEvents_Benchmarks();
  
  PreLoadBegin(PETSC_TRUE, "my preload");
  
  ierr = HornerScheme(1e6, &x); CHKERRQ(ierr);
//  ierr = RandomAccess(1e3); CHKERRQ(ierr);
  ierr = MultipleHornerSchemes(4,1e6); CHKERRQ(ierr);
  
  PreLoadEnd();
  PetscLogStageRegister(&stage, "Benchmarking");
  PetscLogStagePush(stage);
  
//  4294967286/10
  ierr = HornerScheme(1e7, &x); CHKERRQ(ierr);
//  ierr = RandomAccess(1e7); CHKERRQ(ierr);
  ierr = MultipleHornerSchemes(4,1e7); CHKERRQ(ierr);
  
  PetscLogStagePop();
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "RandomAccess"
PetscInt EVENT_RandomAccess;
PetscErrorCode RandomAccess( PetscInt s )
{
  PetscErrorCode ierr;
  Vec v1, v2;
  PetscReal val, *v1p;
  PetscRandom rnd;
  PetscInt *indx;
    
  PetscFunctionBegin;
  ierr = PetscMalloc(s*sizeof(PetscInt), &indx); CHKERRQ(ierr);
  ierr = VecCreateSeq(PETSC_COMM_SELF,s, &v1); CHKERRQ(ierr);
  ierr = VecDuplicate(v1, &v2); CHKERRQ(ierr);
  ierr = PetscRandomCreate(PETSC_COMM_SELF, &rnd); CHKERRQ(ierr);
  ierr = PetscRandomSetType(rnd,PETSCRAND); CHKERRQ(ierr);
  ierr = VecSetRandom(v1, rnd); CHKERRQ(ierr);
  //Time VecCopy
  ierr = VecCopy(v1, v2); CHKERRQ(ierr);
  //Time VecSet
  ierr = VecSet(v2, 0.); CHKERRQ(ierr);
  //Now time random indexing
  ierr = PetscRandomSetInterval(rnd,0,s); CHKERRQ(ierr);
  for( int i = 0; i < s; i++ )
  {
    PetscRandomGetValueReal(rnd, &val);
    indx[i] = (int)val;
  }
  VecGetArray(v1, &v1p);
  PetscLogEventBegin(EVENT_RandomAccess,0,0,0,0);
   
  
  ierr = VecSetValues(v2, s, indx, v1p, INSERT_VALUES); CHKERRQ(ierr);
  
  PetscLogFlops(s);
  PetscLogEventEnd(EVENT_RandomAccess,0,0,0,0);
  
  VecRestoreArray(v1, &v1p);
  ierr = VecDestroy(v1); CHKERRQ(ierr);
  ierr = VecDestroy(v2); CHKERRQ(ierr);
  ierr = PetscRandomDestroy(rnd); CHKERRQ(ierr);
  ierr = PetscFree(indx); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "HornerScheme"
PetscInt EVENT_HornerScheme;
PetscErrorCode HornerScheme(unsigned int iter, double *X)
{
  PetscErrorCode ierr;
  unsigned int i;
  double x = 1.;
  const double me =PETSC_PI * PETSC_MACHINE_EPSILON, 
               x0 = PETSC_PI + me;
  
  PetscFunctionBegin;
  PetscLogEventBegin(EVENT_HornerScheme,0,0,0,0);
  
  for( i = 0; i < iter; i++ )
  {
    x = me + x0 * x;
    x = me + x0 * x;
    x = me + x0 * x;
    x = me + x0 * x;
  }
  
  *X = x;
    
  PetscLogFlops(4 * 2 * iter);
  PetscLogEventEnd(EVENT_HornerScheme,0,0,0,0);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MultipleHornerSchemes"
PetscInt EVENT_MultipleHornerSchemes;
PetscErrorCode MultipleHornerSchemes(int n, unsigned int iter)
{
  PetscErrorCode ierr;
  unsigned int i;
  int j;
  double *x, *me, *x0;
                 
  PetscFunctionBegin;
  ierr = PetscMalloc(n*sizeof(double), &x); CHKERRQ(ierr);
  ierr = PetscMalloc(n*sizeof(double), &x0); CHKERRQ(ierr);
  ierr = PetscMalloc(n*sizeof(double), &me); CHKERRQ(ierr);
  for( j = 0; j < n; j++)
  {
    x[j] = i;
    me[j] = PETSC_PI * PETSC_MACHINE_EPSILON;
    x0[j] = PETSC_PI + me[j]; 
  }
  PetscLogEventBegin(EVENT_MultipleHornerSchemes,0,0,0,0);
  
  for( i = 0; i < iter; i++ )
  {
//    for( j = 0; j < n; j++)
//      x[j] = me[j] + x0[j] * x[j];
    x[0] = me[0] + x0[0] * x[0];
    x[1] = me[1] + x0[1] * x[1];
    x[2] = me[2] + x0[2] * x[2];
    x[3] = me[3] + x0[3] * x[3]; 
  }
  
    
  PetscLogFlops(n * 2 * iter);  
  PetscLogEventEnd(EVENT_MultipleHornerSchemes,0,0,0,0);
  ierr = PetscFree(x); CHKERRQ(ierr);
  ierr = PetscFree(x0); CHKERRQ(ierr);
  ierr = PetscFree(me); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "RegisterEvents_Benchmarks"
PetscErrorCode RegisterEvents_Benchmarks(  )
{
  PetscFunctionBegin;
  PetscLogEventRegister(&EVENT_MultipleHornerSchemes,"MultipleHornerSchemes", 0);
  PetscLogEventRegister(&EVENT_RandomAccess,"RandomAccess", 0);
  PetscLogEventRegister(&EVENT_HornerScheme,"HornerScheme", 0);
  
  PetscFunctionReturn(0);
}

int RunCheck() {return 0;}