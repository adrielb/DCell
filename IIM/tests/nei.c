#include "ImmersedInterfaceMethod.h"

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc, char **args)
{
  PetscErrorCode  ierr;
  PetscFunctionBegin;
  ierr = PetscInitialize(&argc, &args, (char *) 0, ""); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "Start %s\n", __FILE__); CHKERRQ(ierr);

  int d1 = 64;
  Coor d = {1./(d1-1), 1./(d1-1), 1./(d1-1)};
  iCoor size = {d1,d1,0};

  LevelSet ls;
  ierr = LevelSetCreate( size, &ls); CHKERRQ(ierr);
  ierr = GridSetDx(ls->g,d); CHKERRQ(ierr);
  ierr = LevelSetInitializeToStar(ls); CHKERRQ(ierr);
  ierr = VecWrite(ls->g->v, "ls",0); CHKERRQ(ierr);
  ierr = IrregularNodeListWrite(ls->irregularNodes,0); CHKERRQ(ierr);

  PetscReal mu=1;
  IIM iim;
  ierr = IIMCreate(&mu,3,size,20,&iim); CHKERRQ(ierr);

  int len;
  IrregularNode *n, *nodes[iim->Np];

  printf("{");
  for (int j = 0; j < ArrayLength(ls->irregularNodes); ++j) {
    ierr = ArrayGet(ls->irregularNodes,j,(void*)&n); CHKERRQ(ierr);
    ierr = IIMSurfaceNeighbors_2D( iim, ls->g->n,ls->irregularNodeGrid, n, (IrregularNode**)&nodes, &len ); CHKERRQ(ierr);
    if(len>0) printf("{{%f,%f}",nodes[0]->x+1*nodes[0]->ox,nodes[0]->y+1*nodes[0]->oy);
    for (int i = 1; i < len; ++i) {
      printf(",{%f,%f}",nodes[i]->x+nodes[i]->ox,nodes[i]->y+nodes[i]->oy);
    }
    printf("},\n");
  }
  printf("}");

  ierr = LevelSetDestroy(ls); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "End %s\n", __FILE__); CHKERRQ(ierr);
  ierr = PetscFinalize(); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
