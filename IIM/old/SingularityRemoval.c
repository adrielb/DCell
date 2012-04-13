#include "LevelSetMethod.h"

int main(int argc, char **args)
{
  PetscErrorCode  ierr;
  
  ierr = PetscInitialize(&argc, &args, (char *) 0, ""); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "Start\n"); CHKERRQ(ierr);
  
  int d1 = 200, d2 = 200;
  LevelSet2D ls;
  CreateLevelSet2D(d1,d2,&ls);
//  LevelSetInitializeToCircle(ls, PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE);
  LevelSetInitializeToStar(ls,  PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE);
  
  Grid2D rhs;
  CreateGrid2D(d1,d2,&rhs);
  
  IrregularNode *n;
  int ni, nj, k;
  const int nei[5][2] = {1,0, 0,1, -1,0, 0,-1, 0,0};
  PetscReal ox, oy;
  for (int i = 0; i < ls->irregularNodes->len; ++i)
  {
    n = &g_array_index(ls->irregularNodes, IrregularNode, i);
    
    for (int k = 0; k < 5; ++k)
    {
      ni = n->x+nei[k][0];
      nj = n->y+nei[k][1];
      OrthogonalProjection( ls->g2d->v2, ni, nj, &ox, &oy  );
      n->utilde[k] = Bilinear2D( GridFunction2D_Curv, ls->g2d, ni+ox, nj+oy );
      
      /* Correction terms
      C = 0;
      if( -1 * phi[ni][nj] * phi[i][j] )
      {
        C += gamma[k] * PetscSqr( phi[ni][nj] / () ) * (F[ni][nj] - F[i][j]) / 2;
      }
      */
    }
  }

  
  /* print subirregular nodes
  for (int i = 0; i < ls->irregularNodes->len; ++i)
  {
    n = &g_array_index(ls->irregularNodes, IrregularNode, i);
    
    for (k = 0; k < 4; ++k)
    {
      printf("{%d,%d},",n->x+nei[k][0],n->y+nei[k][1]);
    }
  }*/
  
  PetscReal lut; // laplace( u-tilde )
  PetscReal luh; // laplace( u-hat )
    
  for (int i = 0; i < ls->irregularNodes->len; ++i)
  {
    n = &g_array_index(ls->irregularNodes, IrregularNode, i);
    
    // H * laplace( u-tilde )
    lut = 0;
    if( n->sign > 0.)
    {
      for (k = 0; k < 4; ++k)
      {
        lut += n->utilde[k];
      }
      lut -= 4 * n->utilde[4];
    }
    
    // laplace( u-hat ) = laplace( H * u-tilde )
    luh = 0;
    for (k = 0; k < 4; ++k)
    {
      ni = n->x + nei[k][0];
      nj = n->y + nei[k][1];
      if( ls->g2d->v2[ni][nj] > 0 )
      {
        luh += n->utilde[k];
      }
    }
    if( n->sign > 0 )
    {
      luh -= 4 * n->utilde[4];
    }
    
    rhs->v2[n->x][n->y] = -lut + luh; 
    
  }
  
  
  WriteVector("rhs",rhs->v);
  
  WriteVector("phi",ls->g2d->v);
  
  
  
  ierr = PetscPrintf(PETSC_COMM_WORLD, "End\n"); CHKERRQ(ierr);
  ierr = PetscFinalize(); CHKERRQ(ierr);
  return 0; 
}