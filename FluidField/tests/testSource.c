#include "FluidField.h"

void DetermineTimeStep(FluidField f, PetscReal *dt);
PetscErrorCode GrowthTest( FluidField f );

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc, char **args)
{
  PetscErrorCode  ierr;
  ierr = PetscInitialize(&argc, &args, (char *) 0, ""); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "Start %s\n", __FILE__); CHKERRQ(ierr);

                                   int d1 = 32;
  ierr = PetscOptionsSetValue("-da_grid_x","32"); CHKERRQ(ierr);
  ierr = PetscOptionsSetValue("-da_grid_y","32"); CHKERRQ(ierr);
//  ierr = PetscOptionsSetValue("-da_grid_z","32"); CHKERRQ(ierr);
  
  iCoor CELL_CENTER = {0,0,0};
  iCoor s = {d1, d1, 0};
  int lsLen;
  LevelSet ls[2];
  lsLen = 1;
  LevelSetCreate(s, &ls[0]);
  LevelSetInitializeToStar(ls[0]);
  LevelSet new;
  LevelSetCreate(s, &new);
  ierr = VecCopy(ls[0]->g->v, new->g->v); CHKERRQ(ierr);  
  
  PetscReal mu = 1;
  IIM iim;
  IIMCreate( &mu, s, 32, &iim);
  IIMSetForceComponents(iim,ForceComponentNormalSurfaceTension,ForceComponentTangentialSurfaceTension);
//  IIMSetForceComponents(iim,ForceComponentNormalSimple,ForceComponentTangentialSimple);
  
  FluidField f;
  ierr = FluidFieldCreate( &f); CHKERRQ(ierr);
    
  PetscReal ***u, ***v, ***w, dt = 1;
  char file[32];
  sprintf(file, "ls.%d",0);
  WriteVector(file, ls[0]->g->v);
  
  Vec source;
  ierr = FluidFieldGetSource(f, &source); CHKERRQ(ierr);

  for( int t = 1; t < 300; ++t)
  {
    ierr = VecZeroEntries(source); CHKERRQ(ierr);
    for (int j = 0; j < lsLen; ++j)
    {
      iCoor g, m;
      int idx[25];
      PetscReal val[25];
      LevelSetFindMin( ls[j], &m);
      for (int a = 0; a < 5; ++a) {
        for (int b = 0; b < 5; ++b) {
          g.x = m.x + a;
          g.y = m.y + b;
          idx[a+5*b] = LevelSetLocalCoorToGlobalIndex( ls[j], f->da, g );
          val[a+5*b] = -10;
        }
      }
      ierr = VecSetValues(source,25,&idx,&val,ADD_VALUES); CHKERRQ(ierr);
    }  
    ierr = FluidFieldStep(f, iim, lsLen, ls); CHKERRQ(ierr);

    DetermineTimeStep( f, &dt);
    dt = PetscMin(dt, 1);
printf( "%d\t%f\n", t, dt);

    ierr = DAVecGetArray(f->da,f->u,&u); CHKERRQ(ierr);
    ierr = DAVecGetArray(f->da,f->v,&v); CHKERRQ(ierr);
    if( f->is3D )
    {
      ierr = DAVecGetArray(f->da,f->w,&w); CHKERRQ(ierr);
    }
    for (int j = 0; j < lsLen; ++j)
    {
      if( f->is3D )
      {
        ierr = LevelSetAdvectMAC_3D( ls[j], f->d, dt, u, v, w, new ); CHKERRQ(ierr);
      } else {
        ierr = LevelSetAdvectMAC( ls[j], f->d, dt, u, v, new ); CHKERRQ(ierr);
      }
      
      ierr = VecCopy(new->g->v, ls[j]->g->v); CHKERRQ(ierr);
      ierr = IrregularNodeListUpdate(CELL_CENTER, ls[j]); CHKERRQ(ierr);
//      ierr = LevelSetReinitialize( ls[j] ); CHKERRQ(ierr);
      sprintf(file, "ls.%d",t);
      WriteVector(file, ls[j]->g->v);
    }
    ierr = DAVecRestoreArray(f->da,f->u,&u); CHKERRQ(ierr);
    ierr = DAVecRestoreArray(f->da,f->v,&v); CHKERRQ(ierr);
    if( f->is3D )
    {
      ierr = DAVecRestoreArray(f->da,f->w,&w); CHKERRQ(ierr);
    }
  }
  
  ierr = FluidFieldDestroy(f); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "End %s\n", __FILE__); CHKERRQ(ierr);
  ierr = PetscFinalize(); CHKERRQ(ierr);
}

void DetermineTimeStep(FluidField f, PetscReal *dt)
{
  PetscReal max_u, min_u, max_v, min_v, max_w, min_w, U, V, W;
  Coor d = f->d;
  PetscReal CFL = 0.5;
  
  //TODO: try one loop to find min/max 
  VecMax(f->u,PETSC_NULL,&max_u);
  VecMin(f->u,PETSC_NULL,&min_u);
  VecMax(f->v,PETSC_NULL,&max_v);
  VecMin(f->v,PETSC_NULL,&min_v);
  if( f->is3D )
  {
    VecMax(f->w,PETSC_NULL,&max_w);
    VecMin(f->w,PETSC_NULL,&min_w);
  } else {
    max_w = 0;
    min_w = 0;
    d.z = 1;
  }
  
  U = PetscMax( PetscAbs(max_u), PetscAbs(min_u) );
  V = PetscMax( PetscAbs(max_v), PetscAbs(min_v) );
  W = PetscMax( PetscAbs(max_w), PetscAbs(min_w) );

  *dt = CFL / ( U / d.x + V / d.y + W / d.z ); 
}