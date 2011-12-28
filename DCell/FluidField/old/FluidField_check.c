#include "FluidField.h"
#include "MyCheck.h"

PetscErrorCode VecWrite( char *str, Vec v );

START_TEST( viz_StokesFlow )
{
  PetscErrorCode ierr;
  
  
}
END_TEST

START_TEST( viz_FluidFieldStep )
{
  FluidField f;
  DALocalInfo info;
  PetscErrorCode ierr;
  ierr = PetscOptionsSetValue("-da_grid_x","15"); CHKERRQ(ierr);
  ierr = PetscOptionsSetValue("-da_grid_y","15"); CHKERRQ(ierr);
  
  ierr = FluidFieldCreate(&f); CHKERRQ(ierr);
  ierr = FluidFieldStep(f); CHKERRQ(ierr);
  VecWrite("p",f->p);
  VecWrite("px",f->px);
  VecWrite("py",f->py);
  VecWrite("divS",f->div);
  VecWrite("phi",f->phi);
  VecWrite("u",f->u);
  VecWrite("v",f->v);
  
  ierr = DAGetLocalInfo(f->da,&info); CHKERRQ(ierr);
  Divergence(info,f->dx,f->dy,f->u,f->v,f->div);
  VecWrite("divN",f->div);
  ierr = FluidFieldDestroy(f); CHKERRQ(ierr);
}
END_TEST

START_TEST( viz_PoissonNeumanBC )
{
  FluidField f;
  DALocalInfo info;
  PetscReal **div;
  PetscErrorCode ierr;
  ierr = PetscOptionsSetValue("-da_grid_x","15"); CHKERRQ(ierr);
  ierr = PetscOptionsSetValue("-da_grid_y","15"); CHKERRQ(ierr);
  
  ierr = FluidFieldCreate(&f); CHKERRQ(ierr);
  ierr = DAGetLocalInfo(f->da,&info); CHKERRQ(ierr);
  
  ierr = DAVecGetArray(f->da, f->div, &div); CHKERRQ(ierr);
  //div[1*info.my / 4][1*info.mx / 4] = 100;
  //div[2*info.my / 3][2*info.mx / 3] = 100;
  div[1][1] = 100000;
  ierr = DAVecRestoreArray(f->da,f->div,&div); CHKERRQ(ierr);
  PetscReal sum;
  ierr = VecSum(f->div,&sum); CHKERRQ(ierr);
  ierr = VecShift(f->div,-sum/(15*15)); CHKERRQ(ierr);
  
  KSPSetFromOptions(f->kspPhi);
  KSPSolve(f->kspPhi,f->div,f->phi);
  
  VecWrite("phi",f->phi);
  VecWrite("div",f->div);
  ierr = FluidFieldDestroy(f); CHKERRQ(ierr);
}
END_TEST

START_TEST( viz_MomentumDiffusion )
{
  int i,j;
  FluidField f;
  PetscReal q;
  char strU[16];
  DALocalInfo info;
  PetscReal **p;
  PetscErrorCode ierr;
  ierr = PetscOptionsSetValue("-da_grid_x","50"); CHKERRQ(ierr);
  ierr = PetscOptionsSetValue("-da_grid_y","50"); CHKERRQ(ierr);
  
  ierr = FluidFieldCreate(&f); CHKERRQ(ierr);
  ierr = DAGetLocalInfo(f->da,&info); CHKERRQ(ierr);
     q = 1/(f->mu,f->dt);
  ierr = DAVecGetArray(f->da, f->p, &p); CHKERRQ(ierr);
  for (j = 0; j < info.my; ++j)
  {
    for (i = 0; i < info.mx; ++i)
    {
      p[j][i] = 1*i/((double)info.mx);
    }
  }
  ierr = DAVecRestoreArray(f->da,f->p,&p); CHKERRQ(ierr);
  ierr = Gradient2D(info, f->dx,f->dy,f->p,f->px,f->py); CHKERRQ(ierr);
  
  for (int i = 0; i < 100; ++i)
  {
    ierr = VecAYPX(f->u,q,f->px); CHKERRQ(ierr);
    ierr = KSPSolve(f->kspU,f->u,f->u); CHKERRQ(ierr);
    sprintf(strU,"u.%d",i);
    ierr = VecWrite(strU,f->u); CHKERRQ(ierr);    
  }
  
  ierr = FluidFieldDestroy(f); CHKERRQ(ierr);
}
END_TEST

START_TEST( viz_Matrices )
{
  FluidField f;
  PetscErrorCode ierr;
  ierr = PetscOptionsSetValue("-da_grid_x","8"); CHKERRQ(ierr);
  ierr = PetscOptionsSetValue("-da_grid_y","8"); CHKERRQ(ierr);
  
  ierr = FluidFieldCreate(&f); CHKERRQ(ierr);
  
  int s = 256;
  char tmp_dir[256], filename[256];
  PetscViewer view;
  ierr = PetscGetTmp(PETSC_COMM_SELF,tmp_dir,s); CHKERRQ(ierr);
  sprintf(filename,"%s/matU.dat",tmp_dir);
  ierr = PetscViewerASCIIOpen(PETSC_COMM_SELF, filename, &view); CHKERRQ(ierr);
  ierr = MatView(f->matU, view); CHKERRQ(ierr);
  ierr = PetscViewerDestroy(view); CHKERRQ(ierr);
  
  sprintf(filename,"%s/matP.dat",tmp_dir);
  ierr = PetscViewerASCIIOpen(PETSC_COMM_SELF, filename, &view); CHKERRQ(ierr);
  ierr = MatView(f->matP, view); CHKERRQ(ierr);
  ierr = PetscViewerDestroy(view); CHKERRQ(ierr);
  
  sprintf(filename,"%s/matV.dat",tmp_dir);
  ierr = PetscViewerASCIIOpen(PETSC_COMM_SELF, filename, &view); CHKERRQ(ierr);
  ierr = MatView(f->matV, view); CHKERRQ(ierr);
  ierr = PetscViewerDestroy(view); CHKERRQ(ierr);
}
END_TEST

START_TEST( viz_Divergence )
{
  int i,j, X = 70, Y = 70;
  Vec US, VS, DIV;
  PetscReal lx=1, ly=1;
  PetscReal **us, **vs, **div;
  PetscReal xU,yU,xV,yV,dx=lx/(X-2.),dy=ly/(Y-2.),pi=PETSC_PI;
  DA da;
  DALocalInfo info;
  PetscErrorCode ierr;
  ierr = DACreate2d(PETSC_COMM_WORLD,DA_NONPERIODIC,DA_STENCIL_STAR,
          X,Y,PETSC_DECIDE,PETSC_DECIDE,1,1,0,0, &da); CHKERRQ(ierr);
  ierr = DAGetLocalInfo(da, &info); CHKERRQ(ierr);
  ierr = DACreateGlobalVector(da,&DIV); CHKERRQ(ierr);
  ierr = VecDuplicate(DIV,&US); CHKERRQ(ierr);
  ierr = VecDuplicate(DIV,&VS); CHKERRQ(ierr);
  ierr = DAVecGetArray(da,US,&us); CHKERRQ(ierr);
  ierr = DAVecGetArray(da,VS,&vs); CHKERRQ(ierr);
  for( j = 0; j < Y; ++j)
  {
    for( i = 0; i < X; ++i)
    {
      xU = (i-0.5) * dx;
      yU = j * dy;
      xV = i * dx;
      yV = (j-0.5) * dy;
      us[j][i] =  2*cos(pi*yU)*sin(pi*xU)*sin(pi*xU)*sin(pi*yU);
      vs[j][i] = -2*cos(pi*xV)*sin(pi*xV)*sin(pi*yV)*sin(pi*yV);
    }
  }
  ierr = DAVecRestoreArray(da,US,&us); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(da,VS,&vs); CHKERRQ(ierr);
  ierr = Divergence(info,dx,dy,US,VS,DIV); CHKERRQ(ierr);
  ierr = VecWrite("US",US); CHKERRQ(ierr);
  ierr = VecWrite("VS",VS); CHKERRQ(ierr);
  ierr = VecWrite("DIV",DIV); CHKERRQ(ierr);
}
END_TEST

START_TEST( viz_Gradient2D )
{
  int i,j, X = 6, Y = 7;
  Vec P, px, py;
  PetscReal **p;
  DA da;
  DALocalInfo info;
  PetscErrorCode ierr;
  
  ierr = DACreate2d(PETSC_COMM_WORLD,DA_NONPERIODIC,DA_STENCIL_STAR,
        X,Y,PETSC_DECIDE,PETSC_DECIDE,1,1,0,0, &da); CHKERRQ(ierr);
  ierr = DAGetLocalInfo(da, &info); CHKERRQ(ierr);
  ierr = DACreateGlobalVector(da, &P); CHKERRQ(ierr);
  ierr = VecDuplicate(P, &px); CHKERRQ(ierr);
  ierr = VecDuplicate(P, &py); CHKERRQ(ierr);
  ierr = DAVecGetArray(da,P,&p); CHKERRQ(ierr);
  for (j = 0; j < info.my; ++j)
  {
    for (i = 0; i < info.mx; ++i)
    {
      p[j][i] = j;
    }
  }
  ierr = DAVecRestoreArray(da,P,&p); CHKERRQ(ierr);
  mark_point();
  ierr = Gradient2D(info, 1, 1, P, px, py); CHKERRQ(ierr);
  ierr = VecWrite("p",P); CHKERRQ(ierr);
  ierr = VecWrite("px",px); CHKERRQ(ierr);
  ierr = VecWrite("py",py); CHKERRQ(ierr);
}
END_TEST



Suite* FluidField_Suite(void)
{
  Suite *s = suite_create ("Fluid Field Check");
  TCase *tc_core = tcase_create("Core");
  
  tcase_add_test( tc_core,  viz_Gradient2D );
  tcase_add_test( tc_core,  viz_Divergence );
  tcase_add_test( tc_core,  viz_Matrices );
/*
  tcase_add_test( tc_core,  viz_MomentumDiffusion );
  tcase_add_test( tc_core,  viz_PoissonNeumanBC );
  tcase_add_test( tc_core,  viz_FluidFieldStep );
  tcase_add_test( tc_core,  viz_StokesFlow );
  */
//  tcase_add_test( tc_core,  );
  suite_add_tcase( s, tc_core);
  return s;
}