#include "LevelSetMethod.h"

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc, char **args)
{
  PetscErrorCode  ierr;
  PetscFunctionBegin;
  ierr = PetscInitialize(&argc, &args, (char *) 0, ""); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "Start\n"); CHKERRQ(ierr);
  
  LevelSet3D ls;
  ierr = LevelSet3DCreate(100,100,100,&ls); CHKERRQ(ierr);
  ierr = LevelSetInitializeToStar3D(ls); CHKERRQ(ierr);
  
  ierr = IrregularNodeListWrite(ls->irregularNodes,0); CHKERRQ(ierr);
  
//  ierr = ReinitializeLevelSet3D(ls); CHKERRQ(ierr);
  
  PetscViewer viewer;
  ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF,"/home/abergman/Research/DCell/temp/reinit.Real64",FILE_MODE_WRITE,&viewer); CHKERRQ(ierr);
  ierr = VecView(ls->g3d->v,viewer); CHKERRQ(ierr);
  ierr = PetscViewerDestroy(viewer); CHKERRQ(ierr);
//  LevelSet3DDestroy(ls);
  
  ierr = PetscPrintf(PETSC_COMM_WORLD, "End\n"); CHKERRQ(ierr);
  ierr = PetscFinalize(); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}