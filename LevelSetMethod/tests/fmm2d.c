#include "LevelSetMethod.h"
#include "LSM_private.h"

void CheckEikonalProperty( Grid g, Grid temp );

PetscErrorCode SingleFMM()
{
  PetscReal hx = 100;
  Coor dh = {1/hx,1/hx, 0};
  PetscReal radius = 1.0;
  Coor center = (Coor){0,0,0};
  LevelSet ls;
  PetscErrorCode  ierr;

  ierr = LevelSetInitializeToStar2D(dh,center,radius,0.5*radius, 5,&ls); CHKERRQ(ierr);
  ierr = LevelSetSetBandWidth(ls, 30); CHKERRQ(ierr);
  ierr = LevelSetReinitialize(ls); CHKERRQ(ierr);

  printf("Band len: %d\n", ArrayLength(ls->band));

  ierr = GridWrite(ls->phi,0); CHKERRQ(ierr);
  ierr = LevelSetWriteIrregularNodeList(ls,0); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode FMMvsSize()
{
  int i, reps = 10;
  PetscLogDouble t1, t2;
  PetscReal hx;
  PetscReal radius = 1.0;
  Coor center = {0,0,0};
  LevelSet ls;
  PetscLogStage stageLoop;
  PetscErrorCode  ierr;

  ierr = PetscLogStageRegister("dx loop", &stageLoop); CHKERRQ(ierr);
  ierr = PetscLogStagePush(stageLoop); CHKERRQ(ierr);

  for (hx = 32; hx < 1500; ++hx) {
    Coor dh = {1./hx,1./hx,0};
    ierr = LevelSetInitializeToStar2D(dh,center,radius,0.5*radius, 5,&ls); CHKERRQ(ierr);
    ierr = PetscGetTime(&t1); CHKERRQ(ierr);
    for ( i = 0; i < reps; ++i) {
      ierr = LevelSetReinitialize(ls); CHKERRQ(ierr);
    }
    ierr = PetscGetTime(&t2); CHKERRQ(ierr);
    iCoor n = ls->phi->n;
    printf("{%3.1f,%d,%d,%d,%f},\n",hx,n.x,n.x*n.y,ArrayLength(ls->band),(t2-t1)/reps);
    ierr = LevelSetDestroy(ls); CHKERRQ(ierr);
  }
  ierr = PetscLogStagePop(); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode FMMvsParticle()
{
  PetscErrorCode  ierr;
  PetscReal dx = 1/32.;
  PetscReal radius = 1.0;
  Coor center = (Coor){0,0,0};
  Coor dh = {dx,dx,0};

  LevelSet ls;
//  ierr = LevelSetInitializeToCircle(dh,center,radius,&ls); CHKERRQ(ierr);
  ierr = LevelSetInitializeToStar2D(dh,center,radius,0.5*radius, 5,&ls); CHKERRQ(ierr);
  ierr = LevelSetWriteIrregularNodeList(ls,1); CHKERRQ(ierr);
  ierr = GridWrite(ls->phi,1); CHKERRQ(ierr);
  ierr = ArrayWrite(ls->band,"band",1); CHKERRQ(ierr);

  CheckEikonalProperty(ls->phi,ls->tmp);
  ierr = GridSetName(ls->tmp,"tmp"); CHKERRQ(ierr);
  ierr = GridWrite(ls->tmp,0); CHKERRQ(ierr);
  
  ierr = LevelSetInitializeParticles( ls ); CHKERRQ(ierr);
  int i;
  PetscLogStage stageLoop;
  for ( i = 2; i < 3; ++i) {
    printf("i: %d\n", i);
    ierr = LevelSetUpdateIrregularNodeList( ls ); CHKERRQ(ierr);
    ierr = LevelSetReinitialize(ls); CHKERRQ(ierr);
    ierr = ParticleLS_ErrorCorrection2D(ls->pls,ls); CHKERRQ(ierr);
    ierr = LevelSetUpdateIrregularNodeList( ls ); CHKERRQ(ierr);
    ierr = ParticleLS_AdjustRadii(ls->pls,ls); CHKERRQ(ierr);
    ierr = GridWrite(ls->phi,i); CHKERRQ(ierr);
    ierr = LevelSetWriteIrregularNodeList(ls,i); CHKERRQ(ierr);
    if( i == 2 ) {
      ierr = PetscLogStageRegister("FMM loop", &stageLoop); CHKERRQ(ierr);
      ierr = PetscLogStagePush(stageLoop); CHKERRQ(ierr);
    }
  }
  ierr = PetscLogStagePop(); CHKERRQ(ierr);
  return 0;
}

void CheckEikonalProperty( Grid g, Grid temp )
{
  int i, j;
  iCoor p,q;
  PetscReal **phi;
  PetscReal **tmp;
  PetscReal Dxm, Dxp, Dym, Dyp, sign;

  GridGet(g,&phi);
  GridGet(temp,&tmp);
  GridGetBounds(g,&p,&q);
  for( j = p.y + 1; j < q.y - 1; ++j)
  {
    for( i = p.x + 1; i < q.x - 1; ++i)
    {
      sign = phi[j][i] > 0 ? 1 : -1;

      Dxm = sign*(phi[j][i] - phi[j][i-1]);
      Dxm = MAX( Dxm, 0 );
      Dxm = Dxm * Dxm;
      Dxp = sign*(phi[j][i+1] - phi[j][i]);
      Dxp = MIN( Dxp, 0 );
      Dxp = Dxp * Dxp;

      Dym = sign*(phi[j][i] - phi[j-1][i]);
      Dym = MAX( Dym, 0 );
      Dym = Dym * Dym;
      Dyp = sign*(phi[j+1][i] - phi[j][i]);
      Dyp = MIN( Dyp, 0 );
      Dyp = Dyp * Dyp;

      tmp[j][i] = Dxm + Dxp + Dym + Dyp;
        /*
        tmp[j][i] =
          MAX(MAX(phi[j][i]-phi[j-1][i],0.),MAX(phi[j][i]-phi[j+1][i],0.))*MAX(MAX(phi[j][i]-phi[j-1][i],0.),MAX(phi[j][i]-phi[j+1][i],0.)) +
          MAX(MAX(phi[j][i]-phi[j][i-1],0.),MAX(phi[j][i]-phi[j][i+1],0.))*MAX(MAX(phi[j][i]-phi[j][i-1],0.),MAX(phi[j][i]-phi[j][i+1],0.));
         */

    }
  }
}

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc, char **args)
{
  PetscErrorCode  ierr;
  PetscFunctionBegin;
  ierr = DCellInit(); CHKERRQ(ierr);

  ierr = SingleFMM(); CHKERRQ(ierr);

  ierr = DCellFinalize(); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
