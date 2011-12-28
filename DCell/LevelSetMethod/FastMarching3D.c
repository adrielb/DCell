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
  int i,j,k,m;
  PetscReal a,b,c,D;
  PetscReal d1,d2,disc;
  PetscReal ***phi;
  FMMNode node;
  const int numNei = 6;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = GridGet(ls->phi,&phi); CHKERRQ(ierr);
  for( m = 0; m < numNei; ++m)
  {
    i = pos.x + STAR[m][0];
    j = pos.y + STAR[m][1];
    k = pos.z + STAR[m][2];

    // if nei already fixed (behind wavefront) skip it
    if( sign*phi[k][j][i] < ls->bandWidth ) continue;
    // else point is on the wavefront

    a =  0;
    b =  0;
    c = -1;

    AccCoef( sign*phi[k][j][i-1] );
    AccCoef( sign*phi[k][j][i+1] );
    AccCoef( sign*phi[k][j-1][i] );
    AccCoef( sign*phi[k][j+1][i] );
    AccCoef( sign*phi[k-1][j][i] );
    AccCoef( sign*phi[k+1][j][i] );

    // Solve the quadratic equation
    // By quadratic eq:
    disc = PetscSqrtScalar( b*b - 4*a*c );
    d1 = (-b + disc) / ( 2*a );
    d2 = (-b - disc) / ( 2*a );

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

