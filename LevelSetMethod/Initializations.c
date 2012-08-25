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
  pos.x = (int)( (center.x - radius ) / dh.x - BUF );
  pos.y = (int)( (center.y - radius ) / dh.y - BUF );
  pos.z = 0;
  size.x = (int)( 2 * (radius/dh.x + BUF) );
  size.y = (int)( 2 * (radius/dh.y + BUF) );
  size.z = 0;
  ierr = LevelSetCreate(dh, pos, size, &ls); CHKERRQ(ierr);
  ierr = GridGet(ls->phi,&phi); CHKERRQ(ierr);
  ierr = GridGetBounds(ls->phi,&p,&q); CHKERRQ(ierr);
  int est = (int)( PETSC_PI*PetscSqr(radius/dh.x+6) - PETSC_PI*PetscSqr(PetscMax(0,radius/dh.x-6)) );
  ierr = ArraySetSize(ls->band,(int)(1.1*est) ); CHKERRQ(ierr);
  ierr = ArraySetSize(ls->band,0); CHKERRQ(ierr);
  est = (int)( 2*PETSC_PI*radius/dh.x*2.6 ); // Circumference estimate
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

  ierr = LevelSetUpdateIrregularNodeList( ls ); CHKERRQ(ierr);
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
  pos.x = (int)( (center.x - radius ) / dh.x - BUF );
  pos.y = (int)( (center.y - radius ) / dh.y - BUF );
  pos.z = (int)( (center.z - radius ) / dh.z - BUF );
  size.x = (int)( 2 * (radius/dh.x + BUF) );
  size.y = (int)( 2 * (radius/dh.y + BUF) );
  size.z = (int)( 2 * (radius/dh.z + BUF) );

  ierr = LevelSetCreate(dh, pos, size, &ls); CHKERRQ(ierr);
  ierr = GridGet(ls->phi,&phi); CHKERRQ(ierr);
  ierr = GridGetBounds(ls->phi,&p,&q); CHKERRQ(ierr);
  int est = (int)( PETSC_PI*PetscSqr(radius/dh.x+6) - PETSC_PI*PetscSqr(PetscMax(0,radius/dh.x-6)) );
  ierr = ArraySetSize(ls->band,(int)(1.1*est) ); CHKERRQ(ierr);
  ierr = ArraySetSize(ls->band,0); CHKERRQ(ierr);
  est = (int)(2*PETSC_PI*radius/dh.x*2.6); // Circumference estimate
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
          band->z = k;
        } // inside
      } // i
    } // j
  } // k
  
  ierr = LevelSetUpdateIrregularNodeList( ls ); CHKERRQ(ierr);
  ierr = LevelSetReinitialize(ls); CHKERRQ(ierr);

  *lset = ls;

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "LevelSetInitializeFromImage"
PetscErrorCode LevelSetInitializeFromImage( LevelSet ls )
{
  int n, i, j, k;
  const int BUF = 2;
  iCoor p,q;
  iCoor *band;
  PetscReal *phi;
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  ierr = GridGetBounds(ls->phi,&p,&q); CHKERRQ(ierr);
  ierr = GridGet(ls->phi,&phi); CHKERRQ(ierr);
  if( ls->phi->is2D ) {
    PetscReal **phi2D = (PetscReal**)phi;
    for (j = p.y + BUF; j < q.y - BUF; ++j)  {
      for (i = p.x + BUF; i < q.x - BUF; ++i) {
        for (n = 0; n < 4; ++n) {
          if( phi2D[j][i] * phi2D[j + STAR[n][1]][i + STAR[n][0]] <= 0 ) {
            ierr = ArrayAppend(ls->band,&band); CHKERRQ(ierr);
            band->x = i;
            band->y = j;
            break;
          } // if irreg
        } // n
      } // i
    } // j
  } else {
    PetscReal ***phi3D = (PetscReal***)phi;
    for (k = p.z + BUF; k < q.z - BUF; ++k)  {
      for (j = p.y + BUF; j < q.y - BUF; ++j)  {
        for (i = p.x + BUF; i < q.x - BUF; ++i) {
          for (n = 0; n < 6; ++n) {
            if( phi3D[k][j][i] * phi3D[k + STAR[n][2]][j + STAR[n][1]][i + STAR[n][0]] <= 0 ) {
              ierr = ArrayAppend(ls->band,&band); CHKERRQ(ierr);
              band->x = i;
              band->y = j;
              band->z = k;
              break;
            } // if irreg
          } // n
        } // i
      } // j
    } // k
  } // 3D

  ierr = LevelSetUpdateIrregularNodeList(ls); CHKERRQ(ierr);
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
  pos.x = (int)( (center.x - radius ) / dh.x - BUF );
  pos.y = (int)( (center.y - radius ) / dh.y - BUF );
  pos.z = 0;
  size.x = (int)( 2 * (radius/dh.x + BUF) );
  size.y = (int)( 2 * (radius/dh.y + BUF) );
  size.z = 0;
  ierr = LevelSetCreate(dh, pos, size, &ls); CHKERRQ(ierr);
  ierr = GridGet(ls->phi,&phi); CHKERRQ(ierr);
  ierr = GridGetBounds(ls->phi,&p,&q); CHKERRQ(ierr);
  for (j = p.y; j < q.y; ++j) {
    for (i = p.x; i < q.x; ++i) {
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

  ierr = LevelSetUpdateIrregularNodeList( ls ); CHKERRQ(ierr);
  ierr = LevelSetReinitialize(ls); CHKERRQ(ierr);

  *lset = ls;

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "LevelSetInitializeToStar3D"
PetscErrorCode LevelSetInitializeToStar3D( Coor dh, Coor center, PetscReal radius, PetscReal amp, PetscReal numPetals, LevelSet *lset )
{
  int i, j, k;
  const int BUF = 12;
  PetscReal ***phi;
  iCoor *band;
  iCoor p, q;
  iCoor pos,size;
  Coor X;
  PetscReal A;
  PetscReal M;
  LevelSet ls;
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  M = sqrt( 1+2*amp/radius ) * radius;
  pos.x = (int)( (center.x - M) / dh.x - BUF );
  pos.y = (int)( (center.y - M) / dh.y - BUF );
  pos.z = (int)( (center.z - M) / dh.z - BUF );
  size.x = (int)( 2 * (M/dh.x + BUF) );
  size.y = (int)( 2 * (M/dh.y + BUF) );
  size.z = (int)( 2 * (M/dh.z + BUF) );
  ierr = LevelSetCreate(dh, pos, size, &ls); CHKERRQ(ierr);
  ierr = GridGet(ls->phi,&phi); CHKERRQ(ierr);
  ierr = GridGetBounds(ls->phi,&p,&q); CHKERRQ(ierr);
  for (k = p.z; k < q.z; ++k) {
    for (j = p.y; j < q.y; ++j) {
      for (i = p.x; i < q.x; ++i) {
        X.x = i * dh.x - center.x;
        X.y = j * dh.y - center.y;
        X.z = k * dh.z - center.z;
        A = radius + amp + amp * cos( numPetals * atan2(X.y, X.x) );
        A = A / radius;
        phi[k][j][i] = sqrt( (X.x*X.x)/A + (X.y*X.y)/A + (X.z*X.z) ) - radius;
        if( k-p.z > BUF/2 && q.z-k > BUF/2 &&
            j-p.y > BUF/2 && q.y-j > BUF/2 &&
            i-p.x > BUF/2 && q.x-i > BUF/2 ) {
          ierr = ArrayAppend( ls->band, &band); CHKERRQ(ierr);
          band->x = i;
          band->y = j;
          band->z = k;
        }
      }
    }
  }
  
  ierr = LevelSetUpdateIrregularNodeList( ls ); CHKERRQ(ierr);
  ierr = LevelSetReinitialize(ls); CHKERRQ(ierr);

  *lset = ls;

  PetscFunctionReturn(0);
}
