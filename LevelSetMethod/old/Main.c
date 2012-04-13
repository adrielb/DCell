#include "FastMarching.h"
//#include <stdio.h>

#include "../Common/Main.h" 

PetscErrorCode PetscMain() 
{
  PetscErrorCode ierr;
  int WIDTH = 8, HEIGHT = 8;
  
  char mask[8][8]={
    1,1,1,1,1,1,1,1,
    1,1,1,1,1,1,1,1,
    1,1,1,1,1,1,1,1,
    1,1,1,1,0,1,1,1,
    1,1,1,1,1,1,1,1,
    1,1,1,1,1,1,1,1,
    1,1,1,1,1,1,1,1,
    0,1,1,1,1,1,1,1
  };
  
  double phi[8][8];
  
  HEAP *q1 = HeapAlloc(100,&MyDoubleCompMAX);
  HEAP *q2 = HeapAlloc(100,&MyDoubleCompMIN);
  
  /*
  BCFromMask(mask, phi, q1, q2);
  
  SolveEikonal(q1, -1, phi);
  SolveEikonal(q2,  1, phi);
  
  CheckEikonalProperty(phi);
  */
  
  HeapFree(q1);
  HeapFree(q2);
  
  PetscFunctionReturn(0);
}