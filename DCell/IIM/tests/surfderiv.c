#include "petsc.h"
#include "ImmersedInterfaceMethod.h"

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc, char **args)
{
  PetscErrorCode  ierr;
  PetscFunctionBegin;
  ierr = PetscInitialize(&argc, &args, (char *) 0, ""); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "Start\n"); CHKERRQ(ierr);
  
  iCoor s = {64, 64, 0};
  PetscReal mu = 1;
  
  IIM iim;
  ierr = IIMCreate( &mu, 3, s, 20, &iim); CHKERRQ(ierr);
  
  LevelSet ls;
  LevelSetCreate(s,&ls);
  LevelSetInitializeToStar(ls);
  
  iCoor CELL_CENTER = {0,0,0};
  IIMUpdateSurfaceQuantities(iim, CELL_CENTER, ls);
  IrregularNodeListWrite(ls->irregularNodes, 0);
  

  /*
  PetscReal *eta, *xi;    // Local Coordinates in iim->lc
  LocalCoor2DGetVecs(iim->lc, &eta, &xi);
  int len, i, j;
  IrregularNode *n, **nodes, copy;
  GArray *array;
  array = g_array_new( FALSE, FALSE, sizeof(IrregularNode) );
  for( int i = 0; i < ls->irregularNodes->len; i++ )
  {
    n = &g_array_index( ls->irregularNodes, IrregularNode, i);
    ierr = KDTreeRange(iim->kdtree, n, &nodes, &len ); CHKERRQ(ierr);
    for (j = 0; j < len; ++j)
    {
      eta[j] = nodes[j]->x + nodes[j]->ox;
      xi[j]  = nodes[j]->y + nodes[j]->oy;
    }
//    LocalCoor2DSetLength(iim->lc, len);
    LocalCoor2DSolve( iim->lc, n);
    for( j = 0; j < len; j++ )
    {
//      LocalCoor2DToArcLength(iim->lc, n, nodes[j]->nx, nodes[j]->ny, j, &nodes[j]->c);
      PetscMemcpy(&copy,nodes[j],sizeof(IrregularNode));
      copy.nx = eta[j];
      copy.ny = xi[j];
      double n1 = nodes[j]->nx;
      double n2 = nodes[j]->ny;
      copy.df1 = n->nx * n1 + n->ny * n2;
      copy.df2 = n->nx * n2 - n->ny * n1;
      g_array_append_val( array, copy );
    }
    IrregularNodeListWrite(array, i+1);
    if (len != 0 )
      g_array_remove_range(array, 0, len);
  }
  
        */
  LevelSetDestroy(ls);
  IIMDestroy(iim);
  
  ierr = PetscPrintf(PETSC_COMM_WORLD, "End\n"); CHKERRQ(ierr);
  ierr = PetscFinalize(); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
