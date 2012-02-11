#include "LevelSetMethod.h"
/* also initialize to other shapes:
 * LevelSetInitZalesaksDisk
 * Ellipse
 * BinaryImageSegmentation
 */

#undef __FUNCT__
#define __FUNCT__ "LevelSetInitializeToCircle"
PetscErrorCode LevelSetInitializeToCircle( Coor dh, Coor center, PetscReal radius, LevelSet *lset )
{
  int i, j;
  const int BUF = 12;
  PetscReal **phi;
  Coor X;
  iCoor p, q;
  iCoor pos, size;
  iCoor *band;
  LevelSet ls;
  PetscErrorCode ierr = 0;
  
  PetscFunctionBegin;
  pos.x = (center.x - radius ) / dh.x - BUF;
  pos.y = (center.y - radius ) / dh.y - BUF;
  size.x = 2 * (radius/dh.x + BUF);
  size.y = 2 * (radius/dh.y + BUF);
  size.z = 0;
  ierr = LevelSetCreate(dh, pos, size, &ls); CHKERRQ(ierr);
  ierr = GridGet(ls->phi,&phi); CHKERRQ(ierr);
  ierr = GridGetBounds(ls->phi,&p,&q); CHKERRQ(ierr);
  int est = PETSC_PI*PetscSqr(radius/dh.x+6) - PETSC_PI*PetscSqr(PetscMax(0,radius/dh.x-6));
  ierr = ArraySetSize(ls->band,1.1*est); CHKERRQ(ierr);
  ierr = ArraySetSize(ls->band,0); CHKERRQ(ierr);
  est = 2*PETSC_PI*radius/dh.x*2.6; // Circumference estimate
  ierr = ArraySetSize(ls->irregularNodes,est); CHKERRQ(ierr);
  for (j = p.y; j < q.y; ++j) {
    for (i = p.x; i < q.x; ++i) {
      X.x = i * dh.x;
      X.y = j * dh.y;
      phi[j][i] = sqrt( PetscSqr(X.x - center.x) +
                        PetscSqr(X.y - center.y) ) - radius;
      if( -ls->PHI_INF < phi[j][i]/dh.x && phi[j][i]/dh.x < ls->PHI_INF ) {
        ierr = ArrayAppend(ls->band,&band); CHKERRQ(ierr);
        band->x = i;
        band->y = j;
      }
    }
  }

  ierr = LevelSetUpdateIrregularNodeList( ls, ls->phi ); CHKERRQ(ierr);
  ierr = LevelSetReinitialize(ls); CHKERRQ(ierr);

  *lset = ls;

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "LevelSetInitializeToSphere"
PetscErrorCode LevelSetInitializeToSphere( Coor dh, Coor center, PetscReal radius, LevelSet *lset )
{
  int i, j, k;
  const int BUF = 12;
  PetscReal ***phi;
  Coor X;
  iCoor p, q;
  iCoor pos, size;
  iCoor *band;
  LevelSet ls;
  PetscErrorCode ierr = 0;
  
  PetscFunctionBegin;
  pos.x = (center.x - radius ) / dh.x - BUF;
  pos.y = (center.y - radius ) / dh.y - BUF;
  pos.z = (center.z - radius ) / dh.z - BUF;
  size.x = 2 * (radius/dh.x + BUF);
  size.y = 2 * (radius/dh.y + BUF);
  size.z = 2 * (radius/dh.z + BUF);

  ierr = LevelSetCreate(dh, pos, size, &ls); CHKERRQ(ierr);
  ierr = GridGet(ls->phi,&phi); CHKERRQ(ierr);
  ierr = GridGetBounds(ls->phi,&p,&q); CHKERRQ(ierr);
  int est = PETSC_PI*PetscSqr(radius/dh.x+6) - PETSC_PI*PetscSqr(PetscMax(0,radius/dh.x-6));
  ierr = ArraySetSize(ls->band,1.1*est); CHKERRQ(ierr);
  ierr = ArraySetSize(ls->band,0); CHKERRQ(ierr);
  est = 2*PETSC_PI*radius/dh.x*2.6; // Circumference estimate
  ierr = ArraySetSize(ls->irregularNodes,est); CHKERRQ(ierr);
  for (k = p.z; k < q.z; ++k) {
    for (j = p.y; j < q.y; ++j) {
      for (i = p.x; i < q.x; ++i) {
        X.x = i * dh.x;
        X.y = j * dh.y;
        X.z = k * dh.z;
        phi[k][j][i] = sqrt( PetscSqr(X.x - center.x) +
                             PetscSqr(X.y - center.y) +
                             PetscSqr(X.z - center.z) ) - radius;
        if( -ls->PHI_INF < phi[k][j][i]/dh.x && phi[k][j][i]/dh.x < ls->PHI_INF ) {
          ierr = ArrayAppend(ls->band,&band); CHKERRQ(ierr);
          band->x = i;
          band->y = j;
          band->y = k;
        } // inside
      } // i
    } // j
  } // k
  
  ierr = LevelSetUpdateIrregularNodeList( ls, ls->phi ); CHKERRQ(ierr);
  ierr = LevelSetReinitialize(ls); CHKERRQ(ierr);

  *lset = ls;

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "LevelSetInitializeFromImage"
PetscErrorCode LevelSetInitializeFromImage( LevelSet ls )
{
  int i, j;
  const int BUF = 2;
  iCoor p,q;
  iCoor *band;
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  ierr = GridGetBounds(ls->phi,&p,&q); CHKERRQ(ierr);
  for (j = p.y + BUF; j < q.y - BUF; ++j)
  {
    for (i = p.x + BUF; i < q.x - BUF; ++i)
    {
      ierr = ArrayAppend(ls->band,(void*)&band); CHKERRQ(ierr);
      band->x = i;
      band->y = j;
    }
  }

  ierr = LevelSetUpdateIrregularNodeList(ls, ls->phi); CHKERRQ(ierr);
  ierr = LevelSetReinitialize(ls); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "LevelSetInitializeToStar2D"
PetscErrorCode LevelSetInitializeToStar2D( Coor dh, Coor center, PetscReal radius, PetscReal amp, PetscReal numPetals, LevelSet *lset )
{
  int i, j;
  PetscReal **phi;
  iCoor *band;
  iCoor p, q;
  iCoor pos,size;
  Coor X;
  const int BUF = 12;
  PetscReal A;
  LevelSet ls;
  PetscErrorCode ierr = 0;
  
  PetscFunctionBegin;
  pos.x = (center.x - radius ) / dh.x - BUF;
  pos.y = (center.y - radius ) / dh.y - BUF;
  pos.z = 0;
  size.x = 2 * (radius/dh.x + BUF );
  size.y = 2 * (radius/dh.y + BUF );
  size.z = 0;
  ierr = LevelSetCreate(dh, pos, size, &ls); CHKERRQ(ierr);
  ierr = GridGet(ls->phi,&phi); CHKERRQ(ierr);
  ierr = GridGetBounds(ls->phi,&p,&q); CHKERRQ(ierr);
  for (j = p.y; j < q.y; ++j) {
    for (i = p.x; i < q.y; ++i) {
      X.x = i * dh.x - center.x;
      X.y = j * dh.y - center.y;
      A = ( radius + amp + amp * cos( numPetals * atan2(X.y,X.x)) ) / radius;
      phi[j][i] = A * sqrt( X.x*X.x + X.y*X.y ) - radius;
//      if( -thres < phi[j][i]/dh.x && phi[j][i]/dh.x < thres )
      if( j-p.y > BUF/2 && q.y-j > BUF/2 &&
          i-p.x > BUF/2 && q.x-i > BUF/2  )
      {
        ierr = ArrayAppend(ls->band,&band); CHKERRQ(ierr);
        band->x = i;
        band->y = j;
      }
    }
  }

  ierr = LevelSetUpdateIrregularNodeList( ls, ls->phi ); CHKERRQ(ierr);
  ierr = LevelSetReinitialize(ls); CHKERRQ(ierr);

  *lset = ls;

  PetscFunctionReturn(0);
}
/*
#undef __FUNCT__
#define __FUNCT__ "LevelSetInitializeToStar3D"
PetscErrorCode LevelSetInitializeToStar3D( LevelSet ls )
{
  PetscErrorCode ierr;
  int i, j, k;
  PetscReal ***phi = ls->phi->v3;
  iCoor n = ls->phi->n;
  PetscReal I = n.x / 2.-.01,      J = n.y / 2.,      K = n.z / 2.;
//  PetscReal radius = MIN( MIN( I, J), K) / 2 - .1;
  iCoor *band;
  
  PetscFunctionBegin;
    
  for (k = 0; k < n.z; ++k)
  {
    for (j = 0; j < n.y; ++j)
    {
      for (i = 0; i < n.x; ++i)
      {
        phi[k][j][i] = 
          sqrt( 
              ( PetscSqr(i - I) + PetscSqr(j - J) ) *
              (.41 * cos( 8 * atan((j - J)/(i - I))) + 1)
          ) + .5*( PetscSqr(k-K) - K ); // - radius + PetscAbs(k - K);

        ierr = ArrayAppend(ls->band,(void*)&band); CHKERRQ(ierr);
        band->x = i;
        band->y = j;
        band->z = k;
      }
    }
  }
  
  ierr = ArrayDestroy(ls->band); CHKERRQ(ierr);
  ierr = ArrayCreate(sizeof(iCoor), ls->phi->n.x*ls->phi->n.y*ls->bandWidth,&ls->band); CHKERRQ(ierr);
  ierr = LevelSetReinitialize(ls); CHKERRQ(ierr);  
    

  PetscFunctionReturn(0);
}
*/
