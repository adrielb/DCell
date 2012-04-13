#include "LevelSetMethod.h"

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
  ierr = LevelSetCreate( size, d, &ls); CHKERRQ(ierr);
  ierr = GridSetDx(ls->g,d); CHKERRQ(ierr);
  ierr = LevelSetInitializeToBall(ls); CHKERRQ(ierr);

  int i = 0;
  ierr = VecWrite(ls->g->v, "ls",i); CHKERRQ(ierr);
  ierr = IrregularNodeListWrite(ls->irregularNodes,0,i); CHKERRQ(ierr);

  Grid c;
  ierr = GridCreate(size,&c); CHKERRQ(ierr);
  iCoor *b;
  printf("len: %d\n",ArrayLength(ls->band));
  for (int j = 0; j < ArrayLength(ls->band); ++j) {
    ierr = ArrayGet(ls->band,j,(void*)&b); CHKERRQ(ierr);
    c->v2[b->y][b->x] = -1;
    if( PetscAbs(ls->g->v2[b->y][b->x]) < 3 )
    {
      c->v2[b->y][b->x] = GridFunction2D_Curv(ls->g->v2,b->y,b->x,d);
      c->v2[b->y][b->x] = (ls->g->v2[b->y][b->x-1] - 2.*ls->g->v2[b->y][b->x] + ls->g->v2[b->y][b->x+1]) / (d.x*d.x);
    }
  }
  ierr = VecWrite(c->v,"curv",i); CHKERRQ(ierr);

  ierr = LevelSetDestroy(ls); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "End %s\n", __FILE__); CHKERRQ(ierr);
  ierr = PetscFinalize(); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
