#include "LevelSetMethod.h"

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc, char **args)
{
  PetscErrorCode  ierr;
  PetscFunctionBegin;
  ierr = PetscInitialize(&argc, &args, (char *) 0, ""); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "Start\n"); CHKERRQ(ierr);

  iCoor s = {63, 63, 0};
  LevelSet ls, new;
  LevelSetCreate(s, &ls);
  
  LevelSetInitializeToStar(ls,PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE);
  
  IrregularNode *n;
  int count = 0;
  PetscReal x, y, b;
  double i,j, N = 100;
  for( int m = 0; m < ls->irregularNodes->len; ++m)
  {
    n = &g_array_index(ls->irregularNodes,IrregularNode,m);
    /*
     * .---+---.
     * |   |   |
     * +---+---+
     * |   |   |
     * .---+---.
     */
    count = 0;
    for( j = 1; j < N; ++j)
    {
      for( i = 1; i < N; ++i)
      {
        x = i / (N+1) + n->x;
        y = j / (N+1) + n->y;
        b = Bilinear2D( GridFunction2D_Identity, ls->g, x, y);
        if( b > 0. ) 
          count++;
      }
    }
    printf("%f\n",count / N);
  }
  
  ierr = PetscPrintf(PETSC_COMM_WORLD, "End\n"); CHKERRQ(ierr);
  ierr = PetscFinalize(); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}