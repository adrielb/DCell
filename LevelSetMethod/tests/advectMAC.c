#include "LevelSetMethod.h"

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc, char **args)
{
  PetscErrorCode  ierr;
  PetscFunctionBegin;
  ierr = PetscInitialize(&argc, &args, (char *) 0, ""); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "Start\n"); CHKERRQ(ierr);
    
  int d1=1024, d2=d1; 
  LevelSet2D ls, new;
  CreateLevelSet2D(d1,d2,&ls);
  CreateLevelSet2D(d1,d2,&new);
  LevelSetInitializeToStar(ls,PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE);
  ReinitializeLevelSet2D(ls);
  VecCopy(ls->g2d->v,new->g2d->v);
  
  Grid2D u,v;
  CreateGrid2D(d1,d2,&u);
  CreateGrid2D(d1,d2,&v);
  VecSet(u->v,1);
  VecSet(v->v,1);
  
  PetscReal dt = .25;
  Coor dh = {1,1,1};

  WriteVec wv;
  WriteVecCreate(ls->g2d->v,"ls",&wv);
  WriteVecToDisk(wv);
  
  for (int i = 0; i < 120; ++i)
  {
    LevelSet2DAdvectMAC(ls, dh, dt, u->v2, v->v2, new);
    VecCopy(new->g2d->v,ls->g2d->v);
    WriteVecToDisk(wv);
//    UpdateIrregularNodeList(0,0,ls);
//    ReinitializeLevelSet2D(ls);
//    IrregularNodeListWrite(ls->irregularNodes,i);
  }
  
  
  DestroyLevelSet2D(ls);
  
  ierr = PetscPrintf(PETSC_COMM_WORLD, "End\n"); CHKERRQ(ierr);
  ierr = PetscFinalize(); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}