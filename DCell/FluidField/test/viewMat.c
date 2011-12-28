#include "FluidField.h"
#include "FluidField_private.h"

int main(int argc, char **args)
{
  PetscErrorCode  ierr;
  ierr = DCellInit(); CHKERRQ(ierr);

  FluidField f;
  ierr = FluidFieldCreate(PETSC_COMM_WORLD,&f); CHKERRQ(ierr);
  iCoor dims = {6,6,0};
  ierr = FluidFieldSetDims(f, dims); CHKERRQ(ierr);

  Coor dh = {1,1,0};
  iCoor pos = {0,0,0};
  Grid mask;
  ierr = GridCreate(dh,pos,dims,1,&mask); CHKERRQ(ierr);
  ierr = VecSet( mask->v, -1. ); CHKERRQ(ierr);
  PetscReal **mask2D;
  ierr = GridGet(mask, &mask2D); CHKERRQ(ierr);
  mask2D[2][2] = -2;
  ierr = FluidFieldSetMask(f,mask); CHKERRQ(ierr);

  ierr = FluidFieldMatAssemble( f ); CHKERRQ(ierr);
//  ierr = MatWrite(f->mat,"J",0); CHKERRQ(ierr);

  ierr = DCellFinalize(); CHKERRQ(ierr);
  return 0;
}

PetscErrorCode testMixDimGrid()
{

  Coor dh = {1,1,1};
  iCoor pos = {0,0,0};
  iCoor size = {10,10,10};
  int dof = 1;
  Grid g;
  PetscErrorCode  ierr;
  ierr = GridCreate(dh,pos,size,dof,&g); CHKERRQ(ierr);
  PetscReal **mask2D, ***mask3D;
  ierr = GridGet(g, &mask2D); CHKERRQ(ierr);
  ierr = GridGet(g, &mask3D); CHKERRQ(ierr);
  ierr = DCellFinalize(); CHKERRQ(ierr);
  exit(0);
  return 0;
}
