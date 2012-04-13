#include "petsc.h"
#include "ImmersedInterfaceMethod.h"

#define PRINTLINE printf("LINE:%d\n",__LINE__);

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc, char **args)
{
  PetscErrorCode  ierr;
  PetscFunctionBegin;
  ierr = PetscInitialize(&argc, &args, (char *) 0, ""); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "Starting %s\n", __FILE__); CHKERRQ(ierr);
  
  int sq = 32;
  iCoor dim = {sq,sq,sq};
  PetscReal mu = 1;
  IIM iim;
  IIMCreate(&mu, dim, 32, &iim);
  KDTreeSetEps(iim->kdtree,1.8);
  
  LevelSet ls;
  LevelSetCreate(dim,&ls);
  LevelSetInitializeToStar(ls);
//  LevelSetInitializeToBall(ls);

  WriteVector("ls", ls->g->v);
  
  iCoor CELL_CENTER = {0,0,0};
PRINTLINE
  ierr = IIMUpdateSurfaceQuantities( iim, CELL_CENTER, ls ); CHKERRQ(ierr);
PRINTLINE
  ierr = IrregularNodeListWrite(ls->irregularNodes, 0); CHKERRQ(ierr);
//  ierr = IIMUpdateSurfaceDerivatives_2D( iim, ls ); CHKERRQ(ierr);
  
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
      PetscMemcpy(&copy,nodes[j],sizeof(IrregularNode));
      g_array_append_val( array, copy );
    }
    IrregularNodeListWrite(array, i+1);
    if (len != 0 )
      g_array_remove_range(array, 0, len);
  }
  
  LevelSetDestroy(ls);
  IIMDestroy(iim);
  
  ierr = PetscPrintf(PETSC_COMM_WORLD, "End\n"); CHKERRQ(ierr);
  ierr = PetscFinalize(); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}