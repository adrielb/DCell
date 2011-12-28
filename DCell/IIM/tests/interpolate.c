#include "ImmersedInterfaceMethod.h"

void MyFunc( int* b[] );

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc, char **args)
{
  PetscErrorCode  ierr;
  
  ierr = PetscInitialize(&argc, &args, (char *) 0, ""); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "Start %s\n", __FILE__); CHKERRQ(ierr);
  
  iCoor dim = {64, 64, 0};
  
  DA da;
  ierr = DACreate2d(PETSC_COMM_WORLD,DA_NONPERIODIC,DA_STENCIL_BOX,
            dim.x, dim.y, PETSC_DECIDE,PETSC_DECIDE,1,1,0,0, &da); CHKERRQ(ierr);
  
  PetscReal mu = 1;
  IIM iim;
  ierr = IIMCreate(&mu, 3,dim, 12, &iim); CHKERRQ(ierr);
  

  LevelSet ls;
  ierr = LevelSetCreate( dim, &ls); CHKERRQ(ierr);
  ierr = LevelSetInitializeToStar(ls); CHKERRQ(ierr);
  ierr = IrregularNodeListWrite(ls->irregularNodes,0); CHKERRQ(ierr);
  
  Vec vec;
  ierr = DACreateGlobalVector(da,&vec); CHKERRQ(ierr);

  int *ga=0;
  ierr = IIMInterfaceVelocity( iim, da, vec, ga, ls); CHKERRQ(ierr);
  
  ierr = IrregularNodeListWrite(ls->irregularNodes, 0); CHKERRQ(ierr);
  
  ierr = PetscPrintf(PETSC_COMM_WORLD, "End\n"); CHKERRQ(ierr);
  ierr = PetscFinalize(); CHKERRQ(ierr);    
}

void MyFunc( int* b[] )
{
  int i;
  for( i = 0; i < 10; i++ )
  {
    printf("%d\t%d\t%d\n",b[i][0],b[i][1],b[i][2]);
  }
}
