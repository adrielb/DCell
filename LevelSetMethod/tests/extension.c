#include "LevelSetMethod.h"

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc, char **args)
{
  PetscErrorCode  ierr;
  PetscFunctionBegin;
  ierr = PetscInitialize(&argc, &args, (char *) 0, ""); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "Start %s\n", __FILE__); CHKERRQ(ierr);
  
  iCoor size = {2024, 2024, 0};
  LevelSet ls;
  ierr = LevelSetCreate( size, &ls); CHKERRQ(ierr);
  ierr = LevelSetInitializeToStar(ls); CHKERRQ(ierr);
  
  IrregularNode *n;
  for (int i = 0; i < ls->irregularNodes->len; ++i)
  {
    n = IrregularNodeListGet(ls, i);
    double X = n->x + n->op.x,
           Y = n->y + n->op.y;
    
    n->Vn = Bilinear2D( GridFunction2D_Curv, ls->g, X, Y );
  }
  
  IrregularNodeListWrite(ls->irregularNodes,0);
  
  ierr = LevelSetExtendVn( ls ); CHKERRQ(ierr);
  
  WriteVector("Vn", ls->Vn->v);
  
  int i, j;
  PetscReal **phi = ls->g->v2,
            **Vn  = ls->Vn->v2;

  Grid z; // grad(Vn).grad(phi) should equal zero
  ierr = GridCreate( ls->g->n, &z); CHKERRQ(ierr);
  iCoor *band = ArrayGetData( ls->bandNeg ); 
  for (int b = 0; b < ArrayLength(ls->bandNeg); ++b)
  {
    i = band[b].x;
    j = band[b].y;
    
    if( PetscSign( phi[j][i] ) > 0 )
    {
      z->v2[j][i] = 
        PetscMin( phi[j][i+1]-phi[j][i], 0. ) * ( Vn[j][i+1] - Vn[j][i] ) +
        PetscMax( phi[j][i]-phi[j][i-1], 0. ) * ( Vn[j][i] - Vn[j][i-1] ) +
        PetscMin( phi[j+1][i]-phi[j][i], 0. ) * ( Vn[j+1][i] - Vn[j][i] ) +
        PetscMax( phi[j][i]-phi[j-1][i], 0. ) * ( Vn[j][i] - Vn[j-1][i] );
    } else {
      z->v2[j][i] = 
        PetscMin( -phi[j][i+1]+phi[j][i], 0. ) * ( Vn[j][i+1] - Vn[j][i] ) +
        PetscMax( -phi[j][i]+phi[j][i-1], 0. ) * ( Vn[j][i] - Vn[j][i-1] ) +
        PetscMin( -phi[j+1][i]+phi[j][i], 0. ) * ( Vn[j+1][i] - Vn[j][i] ) +
        PetscMax( -phi[j][i]+phi[j-1][i], 0. ) * ( Vn[j][i] - Vn[j-1][i] );
    }
  }
  
  WriteVector("dVn-dPhi", z->v);
  PetscViewer view;
  PetscViewerBinaryOpen(MPI_COMM_SELF,"/tmp/dVn-dPhi.Real64",FILE_MODE_WRITE,&view);
  VecView(z->v,view);
  PetscViewerDestroy(view);
  ierr = GridDestroy(z); CHKERRQ(ierr);
  
  ierr = PetscPrintf(PETSC_COMM_WORLD, "End %s\n", __FILE__); CHKERRQ(ierr);
  ierr = PetscFinalize(); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
