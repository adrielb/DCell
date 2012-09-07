// A second-order fast marching eikonal solver, Ricket and Fomel
#include "LevelSetMethod.h"
#include "LSM_private.h"

inline void FMM_Derivative( const PetscReal bw, const PetscReal *p, PetscReal *a, PetscReal *b, PetscReal *c, PetscReal *D );
inline void FMM_Derivative( const PetscReal bw, const PetscReal *p, PetscReal *a, PetscReal *b, PetscReal *c, PetscReal *D )
{
  const PetscReal D2 = (3.*p[0] - 4.*p[1] + p[2]) / 2.;
  const PetscReal D1 = p[0] - p[ 1];
  PetscReal t;

  if( p[1] < bw && p[2] < bw && D2 > 0 ) {
    // 2nd order Dx
    t = (4*p[1]-p[2])/3;
    *a =  9./4.;
    *b = -9./2.*t;
    *c =  9./4.*t*t;
    *D = D2;
  } else if( p[1] < bw && D1 > 0) {
    // 1st order Dx
    *a = 1.;
    *b = -2.*p[1];
    *c = p[1]*p[1];
    *D = D1;
  } else {
    *a = 0;
    *b = 0;
    *c = 0;
    *D = 0;
  }
}

// max( D+x, D-x, 0)^2
inline void AccumulateCoef( const PetscReal bw, const PetscReal phi[7], PetscReal *a, PetscReal *b, PetscReal *c )
{
  const PetscReal *p = &phi[3];
  const PetscReal D2p = (3*p[0] - 4*p[ 1] + p[ 2]) / 2.;
  const PetscReal D2m = (3*p[0] - 4*p[-1] + p[-2]) / 2.;
  const PetscReal D1p = p[0] - p[ 1];
  const PetscReal D1m = p[0] - p[-1];
  PetscReal Dm, am, bm, cm;
  PetscReal Dp, ap, bp, cp;
  PetscReal t;

  if( 0 <= p[1] && p[1] < bw &&
      0 <= p[2] && p[2] < bw && D2p > 0 ) {
    // 2nd order Dx
    t = (4*p[1]-p[2])/3;
    ap =  9./4.;
    bp = -9./2.*t;
    cp =  9./4.*t*t;
    Dp = D2p;
  } else if( 0 <= p[1] && p[1] < bw && D1p > 0) {
    // 1st order Dx
    ap = 1.;
    bp = -2.*p[1];
    cp = p[1]*p[1];
    Dp = D1p;
  } else {
    ap = 0;
    bp = 0;
    cp = 0;
    Dp = 0;
  }

  if( 0 <= p[-1] && p[-1] < bw &&
      0 <= p[-2] && p[-2] < bw && D2m > 0 ) {
    // 2nd order Dx
    t = (4*p[-1]-p[-2])/3;
    am =  9./4.;
    bm = -9./2.*t;
    cm =  9./4.*t*t;
    Dm = D2m;
  } else if( 0 <= p[-1] && p[-1] < bw && D1m  > 0) {
    // 1st order Dx
    am = 1;
    bm = -2.*p[-1];
    cm = p[-1]*p[-1];
    Dm = D1m;
  } else {
    am = 0;
    bm = 0;
    cm = 0;
    Dm = 0;
  }

  if( Dm > Dp ) {
    *a += am;
    *b += bm;
    *c += cm;
  } else {
    *a += ap;
    *b += bp;
    *c += cp;
  }
}

#undef __FUNCT__
#define __FUNCT__ "FMM_PushNeighbors2D"
PetscErrorCode FMM_PushNeighbors2D( LevelSet ls, MemCache mc, Heap heap, int sign, iCoor pos )
{

  int i,j,k;
  PetscReal a,b,c;
  PetscReal d1,disc;
  PetscReal **phi;
  FMMNode node;
  PetscReal phi_local[7];
  int m, ii, jj;
  const int W = 7; // 3rd order stencil width
  const iCoor min = ls->phi->p;
  const iCoor n = ls->phi->n;
  const iCoor max = {min.x+n.x-1, min.y+n.y-1, 0};
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = GridGet(ls->phi,&phi); CHKERRQ(ierr);
  for( k = 0; k < 4; ++k)
  {
    i = pos.x + STAR[k][0];
    j = pos.y + STAR[k][1];

    // if nei already fixed (behind wavefront) skip it
    if( sign*phi[j][i] < ls->bandWidth ) continue;
    // else point is on the wavefront

    a =  0;
    b =  0;
    c = -1;

    // max(D+x, D-x, 0)^2
    for (m = 0; m < W; ++m) {
      ii = i + m - W/2;
      if( ii < min.x || max.x < ii )
        phi_local[m] = ls->PHI_INF;
      else
        phi_local[m] = sign*phi[j][ii];
    }
    AccumulateCoef( ls->PHI_INF, phi_local, &a, &b, &c);

    // max(D+y, D-y, 0)^2
    for (m = 0; m < W; ++m) {
      jj = j + m - W/2;
      if( jj < min.y || max.y < jj )
        phi_local[m] = ls->PHI_INF;
      else
        phi_local[m] = sign*phi[jj][i];
    }
    AccumulateCoef( ls->PHI_INF, phi_local, &a, &b, &c);

    // Solve the quadratic equation
    disc = PetscSqrtScalar( b*b - 4*a*c );
    d1 = (-b + disc) / ( 2*a );
    if( d1 != d1 ) d1 = -b / (2*a);
    if( d1 != d1 ) {
      printf("\n\n%s: NaN\n", __FUNCT__);
      printf("At: [%d, %d],  sign: %d\n", i,j, sign );
      printf("a: %f\n",a);
      printf("b: %f\n",b);
      printf("c: %f\n",c);
      printf("disc: %f\n",b*b - 4*a*c );
      int J,I;
      for ( J = -2; J <= 2; ++J) {
        for ( I = -2; I <= 2; ++I) {
          printf("%f  ",phi[j+J][i+I]);
        }
        printf("\n");
      }
      ierr = GridSetName(ls->phi,"FMM-err"); CHKERRQ(ierr);
      ierr = GridWrite(ls->phi,0); CHKERRQ(ierr);
      ierr = LevelSetWriteIrregularNodeList(ls,0); CHKERRQ(ierr);
      ierr = LevelSetDestroy(ls); CHKERRQ(ierr);
      exit(1);
      //SETERRQ(0,"FMM NaN");
    } // if NaN

    // add new FMMNode to heap
    ierr = MemCacheAlloc(mc, &node); CHKERRQ(ierr);
    node->phi = d1;
    node->pos.x = i;
    node->pos.y = j;
    ierr = HeapInsert(heap, node); CHKERRQ(ierr);
//    maxHeapSize = maxHeapSize < HeapSize(heap) ? HeapSize(heap) : maxHeapSize;
  } // for k nei of node
  PetscFunctionReturn(0);
}
