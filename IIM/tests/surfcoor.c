#include "ImmersedInterfaceMethod.h"
#include "IIM_private.h"

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc, char **args)
{
  PetscErrorCode  ierr;
  ierr = DCellInit(); CHKERRQ(ierr);

  LocalCoor lc;
  ierr = LocalCoorCreate( &lc ); CHKERRQ(ierr);

  const int len = 5;
  LocalCoorSetLength(lc, len);
  PetscReal *s, *n;
  LocalCoorGetVecs( lc, &s, &n, 0);


  IrregularNode *nodes[5];
  IrregularNode nodestruct[5] = {
      { .X = { -2, 1.2, 0} },
      { .X = { -1, 1.2, 0} },
      { .X = {  0, 1.2, 0} },
      { .X = {  1, 2.2, 0} },
      { .X = {  2, 5.2, 0} }
  };

  int i;
  for ( i = 0; i < len; ++i) {
    nodes[i] = &nodestruct[i];
  }

  Coor dh = {1,1,1};
  IrregularNode node = {
      .nx = -0.707107,
      .ny =  0.707107,
      .nz =  0,
      .X  = {0,1.2,0}
  };
  node.sx =  node.ny;
  node.sy = -node.nx;
  LocalCoorSolve(lc, dh, &node, nodes);

  ierr = ArrayWrite(lc->coor, 0); CHKERRQ(ierr);

  ierr = LocalCoorDestroy(lc); CHKERRQ(ierr);
  ierr = DCellFinalize(); CHKERRQ(ierr);
  return 0;
}
