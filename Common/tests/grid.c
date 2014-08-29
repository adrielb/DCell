#include "Common.h"
#include "Grid.h"
PetscErrorCode PrintGrid(const char title[], Grid g);

int main(int argc, char **args)
{
  PetscErrorCode  ierr;
  ierr = DCellInitialize(&argc,&args,__FILE__); CHKERRQ(ierr);

  Grid g1;
  int dof = 1;
  iCoor gs = {10,11,0};
  Coor dh = {1,1,1};
  iCoor pos = {0,0,0};

  ierr = GridCreate(dh,pos,gs,dof,&g1); CHKERRQ(ierr);
  ierr = VecSet(*((Vec*)g1),2); CHKERRQ(ierr);
  ierr = PrintGrid("GridCreate",g1); CHKERRQ(ierr);
  ierr = GridWrite(g1,0); CHKERRQ(ierr);

  iCoor newpos = {3,4,0};
  ierr = GridResize(g1,newpos,gs); CHKERRQ(ierr);
  ierr = PrintGrid("GridResize shift",g1); CHKERRQ(ierr);
  ierr = GridWrite(g1,1); CHKERRQ(ierr);

  gs.x = 4;
  gs.y = 3;
  ierr = GridResize(g1,newpos,gs); CHKERRQ(ierr);
  ierr = PrintGrid("GridResize size",g1); CHKERRQ(ierr);
  ierr = GridWrite(g1,2); CHKERRQ(ierr);

  gs.x = 8;
  gs.y = 2;
  newpos.x = -3;
  newpos.y = -5;
  ierr = GridResize(g1,newpos,gs); CHKERRQ(ierr);
  ierr = PrintGrid("GridResize shift size",g1); CHKERRQ(ierr);
  ierr = GridWrite(g1,3); CHKERRQ(ierr);

  ierr = GridDestroy(g1); CHKERRQ(ierr);
  ierr = DCellFinalize(); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode PrintGrid(const char title[], Grid g) {
  iCoor p,q;
  int i,j;
  PetscReal **phi;
  PetscErrorCode ierr;

  ierr = GridGet(g,&phi); CHKERRQ(ierr);
  ierr = GridGetBounds(g,&p,&q); CHKERRQ(ierr);
  printf("%s\n   ",title);
  for ( i = p.x; i < q.x; ++i)
    printf("%4d",i);
  printf("\n");
  for ( j = p.y; j < q.y; ++j) {
    printf("%3d:\t",j);
    for ( i = p.x; i < q.x; ++i) {
      printf("%1.1f ", phi[j][i]);
    }
    printf("\n");
  }
  return 0;
}
