#include "LevelSetMethod.h"

PetscErrorCode LevelSetCurvature2D( LevelSet2D ls, Grid2D curvature );

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc, char **args)
{
  PetscErrorCode  ierr;
  PetscFunctionBegin;
  ierr = PetscInitialize(&argc, &args, (char *) 0, ""); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "Start %s\n", __FILE__); CHKERRQ(ierr);
  
  iCoor size = {256,256,64};
  ierr = PetscOptionsGetInt(0,"-width",&size.x,0); CHKERRQ(ierr);
  ierr = PetscOptionsGetInt(0,"-height",&size.y,0); CHKERRQ(ierr);
  ierr = PetscOptionsGetInt(0,"-time",&size.z,0); CHKERRQ(ierr);
  
  Grid imagestack, curvstack;
  ierr = GridCreate(size,&imagestack); CHKERRQ(ierr);
  ierr = GridCreate(size,&curvstack); CHKERRQ(ierr);
  int imageLen = size.x*size.y;
  int stackLen = imageLen*size.z;
  
  int fd;
  int LEN = 256;
  char filename[256];
  ierr = PetscOptionsGetString(0,"-i",filename,LEN,0); CHKERRQ(ierr);
  ierr = PetscBinaryOpen(name,FILE_MODE_READ,&fd); CHKERRQ(ierr);
  ierr = PetscBinaryRead(fd, imagestack->v1, stacklen, PETSC_DOUBLE); CHKERRQ(ierr);
  ierr = PetscBinaryClose(fd); CHKERRQ(ierr);
  
  iCoor lsSize = {size.x, size.y, 0};
  LevelSet ls;
  ierr = LevelSetCreate(lsSize, &ls); CHKERRQ(ierr);
  Grid curv;
  ierr = GridCreate(lsSize, &curv); CHKERRQ(ierr);
  
  for (int i = 0; i < size.z; ++i)
  {
    ierr = PetscMemcpy(ls->g->v1,imagestack->v1[],size.x*size.y); CHKERRQ(ierr);
    ierr = LevelSetCurvature2D(ls,curv); CHKERRQ(ierr);
    ierr = PetscMemcpy(curvstack->v1[],curv->v1,size.x*size.y); CHKERRQ(ierr);
  }
  
  ierr = PetscOptionsGetString(0,"-o",filename,LEN,0); CHKERRQ(ierr);
  ierr = PetscBinaryOpen(name,FILE_MODE_WRITE,&fd); CHKERRQ(ierr);
  ierr = PetscBinaryWrite(fd, curvstack->v1, size.x*size.y*size.z, PETSC_DOUBLE,PETSC_TRUE); CHKERRQ(ierr);
  ierr = PetscBinaryClose(fd); CHKERRQ(ierr);
  
  ierr = LevelSetDestroy(ls); CHKERRQ(ierr);
  ierr = GridDestroy(imagestack); CHKERRQ(ierr);
  ierr = GridDestroy(curvstack); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "End"); CHKERRQ(ierr);
  ierr = PetscFinalize(); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}



#undef __FUNCT__
#define __FUNCT__ "LevelSetCurvature2D"
PetscInt EVENT_LevelSetCurvature2D;
PetscErrorCode LevelSetCurvature2D( LevelSet2D ls, Grid2D curvature )
{
  PetscErrorCode ierr;
  int i, j;
  PetscReal px, py, pxx, pyy, px2, py2, pxy;
  PetscReal **p = ls->g2d->v2;
  PetscReal **k = curvature->v2;
  
  PetscFunctionBegin;
  PetscLogEventBegin(EVENT_LevelSetCurvature2D,0,0,0,0);
  
  for (j = 1; j < ls->g2d->d2 - 1; ++j)
  {
    for( i = 1; i < ls->g2d->d1 - 1; ++i)
    {
      px  = (p[i+1][j]-p[i-1][j]) / two;
      py  = (p[i][j+1]-p[i][j-1]) / two;
      px2 = PetscSqr(px);
      py2 = PetscSqr(py);
      pxx = p[i-1][j] - two*p[i][j] + p[i+1][j];
      pyy = p[i][j-1] - two*p[i][j] + p[i][j+1];
      pxy = (p[i-1][j-1] + p[i+1][j+1] - p[i-1][j+1] - p[i+1][j-1]) / four;
      
      k[i][j]  = pxx*py2 - two*px*py*pxy + pyy*px2;
      k[i][j] /= PETSC_SMALL + sqrt((px2+py2)*(px2+py2)*(px2+py2));
    }
  }
    
  PetscLogEventEnd(EVENT_LevelSetCurvature2D,0,0,0,0);
  PetscFunctionReturn(0);
}