#include "LevelSetMethod.h"

#undef __FUNCT__
#define __FUNCT__ "LevelSet3DCreate"
PetscInt EVENT_LevelSet3DCreate;
PetscErrorCode LevelSet3DCreate( int d1, int d2, int d3, LevelSet3D *ls )
{
  PetscErrorCode ierr;
  LevelSet3D ls3d;
  
  PetscFunctionBegin;
  PetscLogEventBegin(EVENT_LevelSet3DCreate,0,0,0,0);
//  PetscLogEventRegister(&EVENT_LevelSet3DCreate,"LevelSet3DCreate", 0);
  
  PetscNew(struct _LevelSet3D, &ls3d);
  ierr = Grid3DCreate( d1, d2, d3, &ls3d->g3d); CHKERRQ(ierr);
  ls3d->qNeg = HeapAlloc(100,&MyDoubleCompMAX);
  ls3d->qPos = HeapAlloc(100,&MyDoubleCompMIN);  
  ls3d->irregularNodes = g_array_new( FALSE, FALSE, sizeof(IrregularNode) );
  *ls = ls3d;
  
  PetscLogEventEnd(EVENT_LevelSet3DCreate,0,0,0,0);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "LevelSet3DDestroy"
PetscInt EVENT_LevelSet3DDestroy;
PetscErrorCode LevelSet3DDestroy( LevelSet3D ls )
{
  PetscErrorCode ierr;
  
  
  PetscFunctionBegin;
  PetscLogEventBegin(EVENT_LevelSet3DDestroy,0,0,0,0);
//  PetscLogEventRegister(&EVENT_LevelSet3DDestroy,"LevelSet3DDestroy", 0);
  
  ierr = Grid3DDestroy(ls->g3d); CHKERRQ(ierr);
  HeapFree( ls->qNeg );
  HeapFree( ls->qPos);
  g_array_free(ls->irregularNodes,TRUE);
  ierr = PetscFree(ls); CHKERRQ(ierr);
  
  PetscLogEventEnd(EVENT_LevelSet3DDestroy,0,0,0,0);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "LevelSet3DAdvect"
PetscInt EVENT_LevelSet3DAdvect;
PetscErrorCode LevelSet3DAdvect( PetscReal dt, Grid3D gvx, Grid3D gvy, Grid3D gvz, 
    LevelSet3D prev, LevelSet3D new )
{
  PetscErrorCode ierr;
  int i, j, k;
  int d1 = gvx->d1, 
      d2 = gvx->d2,
      d3 = gvx->d3;
  PetscReal px, py, pz;
  PetscReal ***phi = prev->g3d->v3,
            ***fi  =  new->g3d->v3,
            ***vx  = gvx->v3,
            ***vy  = gvy->v3,
            ***vz  = gvz->v3;
  
  PetscFunctionBegin;
  PetscLogEventBegin(EVENT_LevelSet3DAdvect,0,0,0,0);
//  PetscLogEventRegister(&EVENT_LevelSet3DAdvect,"LevelSet3DAdvect", 0);
  
  for (i = 1; i < d1-1; ++i)
  {
    for (j = 1; j < d2-1; ++j)
    {
      for (k = 1; k < d3-1; ++k)
      {
        px = vx[i][j][k] > 0. ? phi[i+1][j][k]-phi[i][j][k] : phi[i][j][k]-phi[i-1][j][k];
        py = vy[i][j][k] > 0. ? phi[i][j+1][k]-phi[i][j][k] : phi[i][j][k]-phi[i][j-1][k];
        pz = vz[i][j][k] > 0. ? phi[i][j][k+1]-phi[i][j][k] : phi[i][j][k]-phi[i][j][k-1];
        fi[i][j][k] = dt * (vx[i][j][k] * px + vy[i][j][k] * py + vz[i][j][k] * pz) + phi[i][j][k];
      }
    }
  }
  printf("bc not implemented\n");
  exit(1);
  PetscLogEventEnd(EVENT_LevelSet3DAdvect,0,0,0,0);
  PetscFunctionReturn(0);
}