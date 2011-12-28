#include "LevelSetMethod.h"

#define a( v ) printf("v:\t%ld\n", ((unsigned long)&n.v) - p)

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc, char **args)
{
  PetscErrorCode  ierr;
  PetscFunctionBegin;
  ierr = DCellInit(); CHKERRQ(ierr);

  int d1 = 64;
  Coor dh = {1./(d1-1), 1./(d1-1), 1./(d1-1)};

  PetscReal radius = 0.5;
  LevelSet ls;
  ierr = LevelSetInitializeToStar2D( dh, (Coor){0,0,0}, radius, 0.1, 7, &ls); CHKERRQ(ierr);
  ierr = GridWrite(ls->phi,0); CHKERRQ(ierr);
  ierr = LevelSetWriteIrregularNodeList(ls->irregularNodes,0); CHKERRQ(ierr);
  /*
  ierr = LevelSetInitializeToCircle( dh, (Coor){0,0,0}, radius, &ls); CHKERRQ(ierr);
  ierr = GridWrite(ls->phi,0); CHKERRQ(ierr);
  ierr = LevelSetWriteIrregularNodeList(ls->irregularNodes,0); CHKERRQ(ierr);
  */
  printf("len: %d\n",ArrayLength(ls->irregularNodes));
  printf("sizeof: %ld\n", sizeof(IrregularNode) );

  IrregularNode n;
  unsigned long p = &n.x;
  a(x);
  a(y);
  a(z);
  printf("\nshift\n");
  a(shift.x);
  a(shift.y);
  a(shift.z);
  printf("\nox\n");
  a(ox);
  a(oy);
  a(oz);
  printf("\n");
  a(nx);
  a(ny);
  a(nz);
  a(sx);
  a(sy);
  a(sz);
  a(rx);
  a(ry);
  a(rz);

  printf("\n");
  a(op.x);
  a(op.y);
  a(op.z);

  IrregularNode *N;
  ierr = ArrayGet(ls->irregularNodes,0,&N); CHKERRQ(ierr);
  printf("[%f, %f, %f] \n",N->ox,N->oy,N->oz);
  printf("[%f, %f, %f] \n",N->nx,N->ny,N->nz);
  printf("[%d, %d, %d] \n",N->x,N->y,N->z);
  printf("[%f, %f, %f] \n",N->op.x,N->op.y,N->op.z);

  double *data = ArrayGetData(ls->irregularNodes);

  for (int i = 0; i < 1e0; ++i) {
    printf("%f,",data[i]);
  }

  ierr = LevelSetDestroy(ls); CHKERRQ(ierr);
  ierr = DCellFinalize(); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
