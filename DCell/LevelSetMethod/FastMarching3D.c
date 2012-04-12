// A second-order fast marching eikonal solver, Ricket and Fomel
#include "LevelSetMethod.h"
#include "LSM_private.h"

// Terrible use of preprocessor directive
#define AccCoef( pp )          \
    D = sign*phi[k][j][i] - pp;\
      if( D > 0 ) {            \
        a += 1;                \
        b += -2*pp;            \
        c += pp*pp;            \
      }

#undef __FUNCT__
#define __FUNCT__ "FMM_PushNeighbors3D"
PetscErrorCode FMM_PushNeighbors3D( LevelSet ls, MemCache mc, Heap heap, int sign, iCoor pos )
{
  int i,j,k,n;
  PetscReal a,b,c;
  PetscReal d1,disc;
  PetscReal ***phi;
  FMMNode node;
  const int numNei = 6;
  const iCoor min = ls->phi->p;
  const iCoor s = ls->phi->n;
  const iCoor max = {min.x+s.x-1, min.y+s.y-1, min.z+s.z-1};
  int m, ii, jj, kk;
  const int W = 7; // 3rd order stencil width
  PetscReal phi_local[W];
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = GridGet(ls->phi,&phi); CHKERRQ(ierr);
  for( n = 0; n < numNei; ++n)
  {
    i = pos.x + STAR[n][0];
    j = pos.y + STAR[n][1];
    k = pos.z + STAR[n][2];

    // if nei already fixed (behind wavefront) skip it
    if( sign*phi[k][j][i] < ls->bandWidth ) continue;
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
        phi_local[m] = sign*phi[k][j][ii];
    }
    AccumulateCoef( ls->PHI_INF, phi_local, &a, &b, &c);

    // max(D+y, D-y, 0)^2
    for (m = 0; m < W; ++m) {
      jj = j + m - W/2;
      if( jj < min.y || max.y < jj )
        phi_local[m] = ls->PHI_INF;
      else
        phi_local[m] = sign*phi[k][jj][i];
    }
    AccumulateCoef( ls->PHI_INF, phi_local, &a, &b, &c);

    // max(D+z, D-z, 0)^2
    for (m = 0; m < W; ++m) {
      kk = k + m - W/2;
      if( kk < min.z || max.z < kk )
        phi_local[m] = ls->PHI_INF;
      else
        phi_local[m] = sign*phi[kk][j][i];
    }
    AccumulateCoef( ls->PHI_INF, phi_local, &a, &b, &c);

    // Solve the quadratic equation
    disc = PetscSqrtScalar( b*b - 4*a*c );
    d1 = (-b + disc) / ( 2*a );
    if( d1 != d1 ) d1 = -b / (2*a);
    if( d1 != d1 ) {
      LINE();
      printf("FMM NaN\n");
//      SETERRQ(0,"FMM NaN");
      exit(1);
    } // if NaN

    // add new FMMNode to heap
    ierr = MemCacheAlloc(mc, &node); CHKERRQ(ierr);
    node->phi = d1;
    node->pos.x = i;
    node->pos.y = j;
    node->pos.z = k;
    ierr = HeapInsert(heap, node); CHKERRQ(ierr);

    // TODO: keep track of max heap size
//    maxHeapSize = maxHeapSize < HeapSize(heap) ? HeapSize(heap) : maxHeapSize;
  } // for m nei of node
  PetscFunctionReturn(0);
}

