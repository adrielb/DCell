#include "ActiveContoursWithoutEdges.h"
#include "LevelSetMethod.h"
#include "FileIO.h"

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
  
  Grid img;
  ierr = GridCreate(size, &img); CHKERRQ(ierr);
  ierr = ReadReal64("-image", img); CHKERRQ(ierr);
  
  ACWE uc;
  ierr = ACWECreate(img, &uc); CHKERRQ(ierr);
  ierr = ReadReal64("-init", uc->ls->g); CHKERRQ(ierr);
  ierr = LevelSetInitializeFromImage(uc->ls); CHKERRQ(ierr);
  ierr = IndexInterior(uc); CHKERRQ(ierr);
  for( int i = 0; i < 100; i++ )
  {
    ierr = SingleStep(uc, img); CHKERRQ(ierr);
    printf("%f,",uc->c1);
    
//    printf("c1: %f\tc2: %f\n", uc->c1,uc->c2);
  }
  ierr = WriteReal64("-o", uc->ls->g); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "End"); CHKERRQ(ierr);
  ierr = PetscFinalize(); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "IndexInterior"
PetscInt EVENT_IndexInterior;
PetscErrorCode IndexInterior(ACWE uc)
{
  int i, c = 0;
  int len = uc->vel->SIZE;
  PetscReal lb, ub;
  PetscReal *ls = uc->ls->g->v1;
  PetscInt *idx;
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  PetscLogEventBegin(EVENT_IndexInterior,0,0,0,0);
  ierr = PetscOptionsGetReal(0,"-lower_bound",&lb,0); CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(0,"-upper_bound",&ub,0); CHKERRQ(ierr);
  
  for( i = 0; i < len; ++i)
  {
    //if( lb < ls[i] && ls[i] < ub ) 
    //{
      c++;
      ArraySetSize(uc->idxArray, c);
      idx = ArrayGetData(uc->idxArray);
      idx[c-1] = i;
    //}
  }  
  
  PetscLogEventEnd(EVENT_IndexInterior,0,0,0,0);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "SingleStep"
PetscInt EVENT_SingleStep;
PetscErrorCode SingleStep( ACWE uc, Grid img )
{
  PetscErrorCode ierr;
  iCoor s= uc->vel->n;
  int x,y;
  PetscInt len = ArrayLength(uc->idxArray), *idx = ArrayGetData(uc->idxArray);
  PetscReal curv, dirac;
  PetscReal *vel = uc->vel->v1,
            *image = img->v1,
            *phi = uc->ls->g->v1,
            maxVel, minVel;
  
  PetscFunctionBegin;
  PetscLogEventBegin(EVENT_SingleStep,0,0,0,0);
  
  ierr = CalculateC1C2(uc->idxArray,uc->ls->g->v1,img->v1, &uc->c1, &uc->c2); CHKERRQ(ierr);
  
  for( int i = 0; i < len; ++i)
  {
    x = idx[i] / s.y;
    y = idx[i] - s.y * x;
    Curvature(uc, x, y, &curv);
    dirac = 1 / (PETSC_PI * ( 1 + phi[idx[i]]*phi[idx[i]])); 
    vel[i] = dirac * ( uc->mu * curv + (
                 uc->l2 * PetscSqr( image[idx[i]] - uc->c2 ) -  
                 uc->l1 * PetscSqr( image[idx[i]] - uc->c1 ) ) ); 
  }
  
  VecMax( uc->vel->v, PETSC_NULL, &maxVel);
  VecMin( uc->vel->v, PETSC_NULL, &minVel);
  maxVel = PetscMax( PetscAbs( maxVel ), PetscAbs( minVel) );
  uc->dt = 10 / maxVel;
  VecAXPY(uc->ls->g->v,uc->dt, uc->vel->v);
  
  PetscLogEventEnd(EVENT_SingleStep,0,0,0,0);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "Curvature"
PetscInt EVENT_Curvature;
PetscErrorCode Curvature(ACWE uc, int x, int y, PetscReal *curv)
{
  int i,j,ii,jj;
  PetscReal p[3][3];
  PetscReal **mask = uc->mask->v2;
  PetscReal **phi = uc->ls->g->v2;
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  PetscLogEventBegin(EVENT_Curvature,0,0,0,0);
//  PetscLogEventRegister(&EVENT_Curvature,"Curvature", 0);
  
  /*
   * It's not dp/dn==0, but eh 
   */
  for( j = -1; j < 2; ++j)
  {
    for( i = -1; i < 2; ++i)
    {
      ii = Clip( x + i );
      jj = Clip( y + j );
      
      if( mask[y+j][x+i] > 0 )
      {
        p[j+1][i+1] = phi[y+j][x+i];
      } else {
        p[j+1][i+1] = phi[y][x];
      }
    }
  }
  
  *curv = GridFunction2D_Curv( (PetscReal**)p,1,1,uc->ls->g->d);
  
  PetscLogEventEnd(EVENT_Curvature,0,0,0,0);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "CalcuateC1C2"
PetscInt EVENT_CalculateC1C2;
PetscErrorCode CalculateC1C2( Array idxArray, PetscReal *phi, PetscReal *image,
  PetscReal *c1, PetscReal *c2)
{
  PetscErrorCode ierr;
  PetscInt *idx = ArrayGetData(idxArray);
  int len = ArrayLength(idxArray);
  const PetscScalar twooverpi = two / PETSC_PI;
  PetscReal  sum1 = zero, sum2 = zero, h;
  
  PetscFunctionBegin;
  PetscLogEventBegin(EVENT_CalculateC1C2,0,0,0,0);
    
  *c1 = *c2 = zero;
  
  for (int i = 0; i < len; ++i)
  {
    h = half * ( one + twooverpi * atan( phi[idx[i]] ) );
    *c1 += image[idx[i]] * h;
    *c2 += image[idx[i]] * (one - h);
    sum1+= h;
    sum2+= (one - h);
  }
  
  *c1 /= sum1;
  *c2 /= sum2;
  
  PetscLogEventEnd(EVENT_CalculateC1C2,0,0,0,0);
  PetscFunctionReturn(0); 
}

#undef __FUNCT__
#define __FUNCT__ "ACWECreate"
PetscErrorCode ACWECreate( Grid img, ACWE *acwe )
{
  PetscErrorCode ierr;
  ACWE uc;
  iCoor s = img->n;
  
  PetscFunctionBegin;
  PetscNew(struct _ACWE, &uc);
  
  uc->dt = one;
  uc->mu = 10 * PetscSqr(1000.);
  uc->l1 = one;
  uc->l2 = one;
  
  ierr = LevelSetCreate(s,&uc->ls); CHKERRQ(ierr);
  ierr = GridCreate(s,&uc->vel); CHKERRQ(ierr);
  ierr = ArrayCreate(sizeof(int),s.x*10,&uc->idxArray); CHKERRQ(ierr);
  *acwe = uc;
  
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ACWEDestroy"
PetscErrorCode ACWEDestroy( ACWE uc )
{
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  ierr = GridDestroy(uc->vel); CHKERRQ(ierr);
  ierr = LevelSetDestroy(uc->ls); CHKERRQ(ierr);
  ierr = PetscFree(uc); CHKERRQ(ierr);
    
  PetscFunctionReturn(0);
}



#undef __FUNCT__
#define __FUNCT__ "RegisterEvents_ACWE"
PetscErrorCode RegisterEvents_ACWE()
{
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  
  PetscFunctionReturn(0);
}