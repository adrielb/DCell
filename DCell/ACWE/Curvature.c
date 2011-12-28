#include "LevelSetMethod.h"
PetscErrorCode Read( Grid g );
PetscErrorCode Write( Grid g );
PetscErrorCode Curvature( LevelSet ls, Grid curv );
inline double GridFunction2D_CurvMag( double **p, int i, int j, Coor d );

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
  Grid curv;
  ierr = GridCreate(size, &curv); CHKERRQ(ierr);
  
  Read( ls->g );
  
  ierr = Curvature(ls, curv); CHKERRQ(ierr);
   
  Write( curv );
  
  ierr = LevelSetDestroy(ls); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "End"); CHKERRQ(ierr);
  ierr = PetscFinalize(); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "Read"
PetscLogEvent EVENT_Read;
PetscErrorCode Read( Grid g )
{
  int fd, LEN = 256;
  char name[256];
  iCoor size = g->n;
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  PetscLogEventBegin(EVENT_Read,0,0,0,0);
//  PetscLogEventRegister("Read", 0, &EVENT_Read);
  ierr = PetscOptionsGetString(0,"-i",name,LEN,0); CHKERRQ(ierr);
  ierr = PetscBinaryOpen(name,FILE_MODE_READ,&fd); CHKERRQ(ierr);
  ierr = PetscBinaryRead(fd, g->v1, size.x*size.y, PETSC_DOUBLE); CHKERRQ(ierr);
  ierr = PetscBinaryClose(fd); CHKERRQ(ierr);
  PetscLogFlops( 0 );
  PetscLogEventEnd(EVENT_Read,0,0,0,0);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "Write"
PetscLogEvent EVENT_Write;
PetscErrorCode Write( Grid g )
{
  int fd, LEN = 256;
  char name[256];
  iCoor size = g->n;
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  PetscLogEventBegin(EVENT_Write,0,0,0,0);
//  PetscLogEventRegister("Write", 0, &EVENT_Write);
  ierr = PetscOptionsGetString(0,"-o",name,LEN,0); CHKERRQ(ierr);
  ierr = PetscBinaryOpen(name,FILE_MODE_WRITE,&fd); CHKERRQ(ierr);
  ierr = PetscBinaryWrite(fd, g->v1, size.x*size.y, PETSC_DOUBLE,PETSC_TRUE); CHKERRQ(ierr);
  ierr = PetscBinaryClose(fd); CHKERRQ(ierr);
  PetscLogFlops( 0 );
  PetscLogEventEnd(EVENT_Write,0,0,0,0);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "Curvature"
PetscLogEvent EVENT_Curvature;
PetscErrorCode Curvature( LevelSet ls, Grid curv )
{
  int i,j;
  iCoor s = ls->g->n;
  double **phi = ls->g->v2, **k = curv->v2;
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  PetscLogEventBegin(EVENT_Curvature,0,0,0,0);
//  PetscLogEventRegister("Curvature", 0, &EVENT_Curvature);
  
  for (j = 1; j < s.y-1; ++j)
  {
    for (i = 1; i < s.x-1; ++i)
    {
      if ( i < 5 && j < 5 ) printf("%2.2f ", phi[j][i]); 
      k[j][i] = GridFunction2D_Curv( phi, i, j, ls->g->d );
      if( k[j][i] != k[j][i] ) k[j][i] = 0;

    }
  }
  
  PetscLogFlops( 0 );
  PetscLogEventEnd(EVENT_Curvature,0,0,0,0);
  PetscFunctionReturn(0);
}

inline double GridFunction2D_CurvMag( double **p, int i, int j, Coor d )
{
  double px, py, px2, py2, pxx, pyy, pxy, mag, k;
  px = GridFunction2D_DerivX( p, i, j, d );
  py = GridFunction2D_DerivY( p, i, j, d );
  px2 = px * px;
  py2 = py * py;
  pxx = (p[j][i-1] - 2.*p[j][i] + p[j][i+1]) / (d.x*d.x);
  pyy = (p[j-1][i] - 2.*p[j][i] + p[j+1][i]) / (d.y*d.y);
  pxy = (p[j-1][i-1] + p[j+1][i+1] - p[j+1][i-1] - p[j-1][i+1]) / (4.*d.x*d.y);
  k   = pxx*py2 - 2.*px*py*pxy + pyy*px2;
  mag = sqrt((px2+py2)*(px2+py2)*(px2+py2));
  k  /= mag;
  return mag;
}