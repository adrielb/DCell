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
  int lsLen=1;
  LevelSet ls[lsLen];
  
  LevelSetCreate(s, &ls[0]);
  LevelSetInitializeToStar(ls[0]);
  
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
  
  int t, j;
  for( t = 1; t < 2; ++t)
  {
    ierr = FluidFieldStep(f, iim, lsLen, ls); CHKERRQ(ierr);

    DetermineTimeStep( f, &dt);
    dt = PetscMin(dt, 1);
    
    printf( "%d\t%f\n", t, dt);

    for (j = 0; j < lsLen; ++j)
    {
      ierr = IIMInterfaceVelocity( iim, f->da, &f->u, &f->gaU, ls[j] ); CHKERRQ(ierr); 
      ierr = LevelSetExtendVn( ls[j] ); CHKERRQ(ierr);
      ierr = LevelSetAdvectVn( ls[j] ); CHKERRQ(ierr);
      
      ierr = IrregularNodeListUpdate(CELL_CENTER, ls[j]); CHKERRQ(ierr);
//      ierr = LevelSetReinitialize( ls[j] ); CHKERRQ(ierr);
      sprintf(file, "ls.%d",t);
      WriteVector(file, ls[j]->g->v);
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