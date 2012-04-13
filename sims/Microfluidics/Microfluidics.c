#include "DWorld.h"

PetscErrorCode SetPressureBC( FluidField f );

int main(int argc, char **args) {
  PetscErrorCode ierr;
  ierr = DCellInit(); CHKERRQ(ierr);

  PetscReal dx = 1;
  iCoor size = {1625,1145,0};
//  iCoor size = {253,341,0};
  int fd;
  Coor dh = {dx,dx,dx};
  iCoor pos = {0,0,0};
  Grid chip;
  ierr = GridCreate(dh,pos,size,1,&chip); CHKERRQ(ierr);
  ierr = PetscInfo(0,"Reading image file\n"); CHKERRQ(ierr);
  ierr = PetscBinaryOpen("/scratch/n/BL Big",FILE_MODE_READ,&fd); CHKERRQ(ierr);
  ierr = PetscBinaryRead(fd,chip->v1,size.x*size.y,PETSC_DOUBLE); CHKERRQ(ierr);
  ierr = PetscBinaryClose(fd); CHKERRQ(ierr);
  ierr = PetscInfo(0,"Writing image file\n"); CHKERRQ(ierr);
  ierr = GridWrite(chip,0); CHKERRQ(ierr);

  FluidField fluid;
  ierr = FluidFieldCreate(PETSC_COMM_WORLD, &fluid);  CHKERRQ(ierr);
  ierr = FluidFieldSetDims(fluid,size); CHKERRQ(ierr);
  ierr = FluidFieldSetDx(fluid,dx); CHKERRQ(ierr);
  ierr = FluidFieldSetMask(fluid, chip); CHKERRQ(ierr);
  ierr = FluidFieldSetup(fluid); CHKERRQ(ierr);
  ierr = SetPressureBC(fluid); CHKERRQ(ierr);
  ierr = KSPSolve(fluid->ksp,fluid->rhs,fluid->vel); CHKERRQ(ierr);
  ierr = FluidFieldWrite( fluid,0); CHKERRQ(ierr);

  ierr = FluidFieldDestroy(fluid); CHKERRQ(ierr);
  ierr = GridDestroy(chip); CHKERRQ(ierr);
  ierr = DCellFinalize(); CHKERRQ(ierr);
  return 0;
}

PetscErrorCode SetPressureBC( FluidField f )
{
  int i,j;
  PetscReal **mask, m;
  PetscReal ***rhs,p;
  DALocalInfo g;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscInfo(0,"Entering SetPressureBC()"); CHKERRQ(ierr);
  ierr = DAGetLocalInfo(f->daB, &g); CHKERRQ(ierr);
  ierr = GridGet(f->mask, &mask); CHKERRQ(ierr);
  ierr = DAVecGetArrayDOF(f->daV,f->rhs,&rhs); CHKERRQ(ierr);
  for (j = g.ys; j < g.ys+g.ym; ++j) {
    for (i = g.xs; i < g.xs+g.xm; ++i) {
      m = mask[j][i];
      if( m > -1.5 ) continue;
      if( PetscAbs(m+2) < .01 ) p = 20;
      if( PetscAbs(m+3) < .01 ) p = 19;
      if( PetscAbs(m+4) < .01 ) p = 10;
      rhs[j][i][CELL_CENTER] = p;
    }
  }
  ierr = DAVecRestoreArrayDOF(f->daV,f->rhs,&rhs); CHKERRQ(ierr);
  ierr = PetscInfo(0,"Exiting SetPressureBC()"); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode CreateChip()
{
  PetscErrorCode ierr;
  PetscReal dx = 1;
  Coor len = {32,16,0};
  Coor dh = {dx,dx,0};
  iCoor size = {len.x/dx,len.y/dx,0};
  printf("MX = %d;\n", size.x);
  printf("MY = %d;\n", size.y);

  iCoor pos = {0,0,0};
  Grid chip;
  ierr = GridCreate(dh,pos,size,1,&chip); CHKERRQ(ierr);
  ierr = VecSet(chip->v,-1); CHKERRQ(ierr);
  int borderwidth = 3;
  ierr = GridDrawBorder( chip, borderwidth, 1 ); CHKERRQ(ierr);
  PetscReal  sq = 2;
  {
    Coor loc = {8,8,0};
    Coor lo = {loc.x-sq,loc.y-sq,0};
    Coor hi = {loc.x+sq,loc.y+sq,0};
    ierr = GridFillRectangle(chip, lo, hi, -2); CHKERRQ(ierr);
  }
  {
    Coor loc = {24,8,0};
    Coor lo = {loc.x-sq,loc.y-sq,0};
    Coor hi = {loc.x+sq,loc.y+sq,0};
    ierr = GridFillRectangle(chip, lo, hi, -3); CHKERRQ(ierr);
  }
  ierr = GridSetName(chip,"chip"); CHKERRQ(ierr);
  ierr = GridWrite(chip,0); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
