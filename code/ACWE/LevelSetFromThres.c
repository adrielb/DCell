#include "LevelSetMethod.h"

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc, char **args)
{
  PetscErrorCode  ierr;
  PetscFunctionBegin;
  ierr = PetscInitialize(&argc, &args, (char *) 0, ""); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "Start %s\n", __FILE__); CHKERRQ(ierr);
  
  iCoor size = {256,256,0};
  ierr = PetscOptionsGetInt(0,"-width",&size.x,0); CHKERRQ(ierr);
  ierr = PetscOptionsGetInt(0,"-height",&size.y,0); CHKERRQ(ierr);
  
  
  LevelSet ls;
  ierr = LevelSetCreate(size, &ls); CHKERRQ(ierr);
  
  int fd;
  int LEN = 256;
  char name[256];
  ierr = PetscOptionsGetString(0,"-i",name,LEN,0); CHKERRQ(ierr);
  ierr = PetscBinaryOpen(name,FILE_MODE_READ,&fd); CHKERRQ(ierr);
  ierr = PetscBinaryRead(fd, ls->g->v1, size.x*size.y, PETSC_DOUBLE); CHKERRQ(ierr);
  ierr = PetscBinaryClose(fd); CHKERRQ(ierr);
    
  ierr = LevelSetInitializeFromImage(ls); CHKERRQ(ierr);
      
  ierr = PetscOptionsGetString(0,"-o",name,LEN,0); CHKERRQ(ierr);
  ierr = PetscBinaryOpen(name,FILE_MODE_WRITE,&fd); CHKERRQ(ierr);
  ierr = PetscBinaryWrite(fd, ls->g->v1, size.x*size.y, PETSC_DOUBLE,PETSC_TRUE); CHKERRQ(ierr);
  ierr = PetscBinaryClose(fd); CHKERRQ(ierr);
  
  ierr = LevelSetDestroy(ls); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "End"); CHKERRQ(ierr);
  ierr = PetscFinalize(); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
