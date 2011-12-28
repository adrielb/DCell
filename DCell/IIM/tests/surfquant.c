#include "ImmersedInterfaceMethod.h"

void InterfacialForceAdhesion( IrregularNode *n );

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc, char **args)
{
  PetscErrorCode  ierr;
  PetscFunctionBegin;
  ierr = PetscInitialize(&argc, &args, (char *) 0, ""); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "Start\n"); CHKERRQ(ierr);

  int d1 = 64;
  Coor d = {1./(d1-1), 1./(d1-1), 1./(d1-1)};
  iCoor size = {d1,d1,0};

  LevelSet ls;
  ierr = LevelSetCreate( size, &ls); CHKERRQ(ierr);
  ierr = GridSetDx(ls->g,d); CHKERRQ(ierr);
  ierr = LevelSetInitializeToBall(ls); CHKERRQ(ierr);
  ierr = VecWrite(ls->g->v, "ls",0); CHKERRQ(ierr);

  PetscReal mu=1;
  IIM iim;
  ierr = IIMCreate(&mu,3,size,20,&iim); CHKERRQ(ierr);
  IIMSetForceComponents(iim,InterfacialForceAdhesion);

  iCoor CELL_CENTER={0,0,0};
  ierr = IIMUpdateSurfaceQuantities(iim,CELL_CENTER,ls); CHKERRQ(ierr);
  ierr = IrregularNodeListWrite(ls->irregularNodes,0); CHKERRQ(ierr);

  ierr = LevelSetDestroy(ls); CHKERRQ(ierr);
  ierr = IIMDestroy(iim); CHKERRQ(ierr);
  
  ierr = PetscPrintf(PETSC_COMM_WORLD, "End\n"); CHKERRQ(ierr);
  ierr = PetscFinalize(); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

void InterfacialForceAdhesion( IrregularNode *n )
{
  const int gx = 0, gy = -1;
  n->f1 = n->nx * gx + n->ny * gy;
  n->f2 = -n->ny * gx + n->nx * gy;
}
