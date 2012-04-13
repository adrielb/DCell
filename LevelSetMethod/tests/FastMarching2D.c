#include "LevelSetMethod.h"
#include "stdlib.h"

#undef __FUNCT__
#define __FUNCT__ "LevelSetReinitialize_2D"
PetscErrorCode LevelSetReinitialize_2D( LevelSet ls, MemCache mc, Heap heap )
{
  double **phi;
  PetscReal dist;
  int len = ArrayLength(ls->irregularNodes);
  IrregularNode *nodes = ArrayGetData(ls->irregularNodes);
  IrregularNode *n;
  PetscErrorCode ierr;  
  
  PetscFunctionBegin;
  ierr = GridGet(ls->phi,&phi); CHKERRQ(ierr);
  /*Prevent adding bc node more than once to heap.
     *(irregnodelist may contain multiple ortho-proj for the same node if list generated from grid intersections)
     */
  for( int i = 0; i < len; i++ ) {
    n = &nodes[i];
    if( PetscAbs( n->sign * phi[n->y][n->x] - ls->PHI_INF ) > 0.1 ) continue;
    dist = n->sign * sqrt( PetscSqr(n->op.x) + PetscSqr(n->op.y) );
    phi[n->y][n->x] = dist;
    if( dist <= 0. )
    {
      ierr = MemCacheAlloc(mc,&node); CHKERRQ(ierr);
      node->
      ierr = HeapInsert(h,node); CHKERRQ(ierr);
      FMMAddNodeToHeap( ls->qNeg, (iCoor){n->x, n->y, 0}, dist);
    }
    if( dist >= 0. )
    {
      FMMAddNodeToHeap( ls->qPos, n->x, n->y, 0, dist);
    }
  }
  

  ierr = FMMSolveEikonal_2D(ls->qPos, 1, ls->bandWidth, ls->band, phi, ls->phi->n, ls->phi->p, ls->PHI_INF); CHKERRQ(ierr);
  ierr = FMMSolveEikonal_2D(ls->qNeg,-1, ls->bandWidth, ls->band, phi, ls->phi->n, ls->phi->p, ls->PHI_INF); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

double CheckLimits(double **phi, PetscReal PHI_INF, iCoor lo, iCoor hi, int i, int j )
{
  if( i == lo.x || i == hi.x ||
      j == lo.y || j == hi.y )
    return PHI_INF;
  return phi[j][i];
}

#undef __FUNCT__
#define __FUNCT__ "FMMSolveEikonal_2D"
PetscErrorCode FMMSolveEikonal_2D(HEAP *q, int sign, PetscReal bandWidth, Array aBand, PetscReal **phi, iCoor s, iCoor p, PetscReal PHI_INF)
{
  int i,k;
  iCoor lo = {p.x-1,   p.y-1,  0};
  iCoor hi = {p.x+s.x, p.y+s.y,0};
  int nei[4][2] = {{1,0},{-1,0},{0,1},{0,-1}};
  const PetscReal TOL = .001; //tolerance when comparing to INF
  iCoor *band;
  int ni, nj;
  double x, y, d;
  FMMNode n;
  PetscErrorCode ierr=0;

  PetscFunctionBegin;
  
i = 0;
  while( HeapSize(q) != 0 )
  {
    n = (FMMNode)(HeapNext(q).v);
        
    //Add node to band
    ierr = ArrayAppend(aBand,(void*)&band); CHKERRQ(ierr);
    band->x = n->i;
    band->y = n->j;
    
    for( k = 0; k < 4; ++k)
    {
      ni = n->i + nei[k][0];
      nj = n->j + nei[k][1];
      if(ni == lo.x || nj == lo.y || ni == hi.x || nj == hi.y ) continue;
      if( ABS( sign*phi[ nj ][ ni ] - PHI_INF) < TOL )  // phi == +/-INF
      {
        x = MIN( ABS( CheckLimits( phi, PHI_INF, lo, hi, ni + 1 , nj ) ),
                 ABS( CheckLimits( phi, PHI_INF, lo, hi, ni - 1 , nj ) ) );
      
        y = MIN( ABS( CheckLimits( phi, PHI_INF, lo, hi, ni, nj + 1 ) ),
                 ABS( CheckLimits( phi, PHI_INF, lo, hi, ni, nj - 1 ) ) );

        if( ABS(x - PHI_INF) < TOL ) // x == +PHI_INF
        {
          d = y + 1;
        } else if( ABS(y - PHI_INF) < TOL ) { // y == +PHI_INF
          d = x + 1;
        } else {
          d = FMMDist2D( x, y);
//          d = Distance( x, y);
          if( x == x && y == y)
          {
            if( d != d ) // Some type of error, debug output:
            {
              printf("\n**************\n");
              printf("(%d,%d)\n",ni,nj);
              printf("%f,%f; %f\n", x, y, d);
              printf("{%f,%f, %f,%f}\n", phi[nj][ni+1], phi[nj][ni-1], phi[nj+1][ni],phi[nj-1][ni]);
              printf("FMM2D, NaN distance");
//              ierr = VecWrite(g->v,"phi.debug",0); CHKERRQ(ierr);
              exit(1);
            }
          }
        }
        phi[nj][ni] = sign * d;
/* DEBUG: output fmm evaluation after every iteration
char w[8];
sprintf(w,"%s.%d","fmm",sign);
ierr = VecWrite(g->v,w,i++); CHKERRQ(ierr);
*/
        if( PetscAbs(n->dist) > bandWidth ) break; // Stop marching if node beyond bandwidth
        FMMAddNodeToHeap(q, ni, nj, 0, phi[nj][ni]);
      }
    }
    free( n ); //TODO: eliminate this free and malloc a block of FMMNodes
  }
  while( HeapSize(q) != 0 ) // Need to free nodes still in heap if dist above narrow band width
  {
    n = (FMMNode)(HeapNext(q).v);
    free( n ); 
  }

  PetscFunctionReturn(0);
}

inline double Distance(double x, double y) 
{
  return 1/2. * ( x + y + sqrt( 2. - (x-y) * (x-y) ) );
}

double FMMDist2D(double x, double y)
{
  return 0.5 * (x + y + sqrt(2. - x*x + 2.*x*y - y*y));
}
/*
void BCFromMask( int WIDTH, int HEIGHT, double **mask, double **phi, HEAP *q1, HEAP *q2 )
{
  printf("BCFromMask DEPRECATED");
  exit(1);
  int i,j,k;
  int nei[4][2] = {{1,0},{-1,0},{0,1},{0,-1}};
  const double dist[5] = {PHI_INF, 0.5, 0.5/sqrt(2), 0.5/sqrt(2), 0.5/sqrt(2)};
  int count = 0;
  int nj, ni;
  double sgn = 1;
  
  HEAP *q;
    
  for( i = 0; i < WIDTH; ++i)
    for( j = 0; j < HEIGHT; ++j)
      phi[i][j] = INF;
  
  for( i = 0; i < WIDTH; ++i )
  {
    for( j = 0; j < HEIGHT; ++j)
    {
      count = 0;
      if( mask[i][j] < 0 )
      {
        sgn = -1;
        q   = q1; 
      } else {
        sgn = 1;
        q   = q2;
      }
      for( k = 0; k < 4; ++k)
      {
        ni = i+nei[k][0];
        nj = j+nei[k][1];
        if( ni == -1 || nj == -1 || ni == WIDTH || nj == HEIGHT ) continue;
        if( sgn * mask[ ni ][ nj ] < 0 )
        {
          ++count;
        }
      }
      if( count > 0 )
      {
        phi[i][j] = sgn * dist[count];
        FMMAddNodeToHeap(q, i, j, 0, phi[i][j]);
      }
    }
  }
//  PrintPhi( WIDTH, HEIGHT, phi );
}
*/
