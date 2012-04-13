#include "FluidField.h"

int main(int argc, char **args)
{
  PetscErrorCode  ierr;
  ierr = PetscInitialize(&argc, &args, (char *) 0, ""); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "Start\n"); CHKERRQ(ierr);
  
                                   int d1 = 5;
  ierr = PetscOptionsSetValue("-da_grid_x","5"); CHKERRQ(ierr);
  
                                   int d2 = 7;
  ierr = PetscOptionsSetValue("-da_grid_y","7"); CHKERRQ(ierr);
  
  FluidField f;
  ierr = FluidFieldCreate(&f); CHKERRQ(ierr);  
  
  iCoor s = {d1,d2,0};
  Grid g;
  ierr = GridCreate(s,&g); CHKERRQ(ierr);
  
  PetscReal *p;
  VecGetArray(f->p, &p);
  
  for (int i = 0; i < g->n.x*g->n.y; ++i)
  {
    g->v1[i] = i;
    p[i] = i;
  }
  
  PetscReal **pp;
  DAVecGetArray(f->da, f->p, &pp);
  for (int j = 0; j < d2; ++j)
  {
    for (int i = 0; i < d1; ++i)
    {
      printf( "%1.0f, %1.0f\t", g->v2[j][i], pp[j][i]);
    }
    printf("\n");
  }
  DAVecRestoreArray(f->da, f->p, &pp);
  ierr = GridDestroy(g); CHKERRQ(ierr);
  ierr = FluidFieldDestroy(f); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "End\n"); CHKERRQ(ierr);
  ierr = PetscFinalize(); CHKERRQ(ierr);
  return 0;
}