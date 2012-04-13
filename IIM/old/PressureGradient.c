#include "ImmersedInterfaceMethod.h"

PetscReal IIMxGradient(IrregularNode *n, PetscReal **p );
PetscReal IIMyGradient(IrregularNode *n, PetscReal **P );

void PressureGradient( )
{
  
  for( k = 1; k < d3; ++k)
  {
    for( j = 1; j < d2; ++j)
    {
      for( i = 1; i < d1; ++i)
      {
        px[k][j][i] = p[k][j][i+1] - p[k][j][i-1];
        py[k][j][i] = p[k][j+1][i] - p[k][j-1][i];
        pz[k][j][i] = p[k+1][j][i] - p[k-1][j][i];
      }
    }
  }
}

#undef __FUNCT__
#define __FUNCT__ "IIMPressureGradient"
PetscInt EVENT_IIMPressureGradient;
PetscErrorCode IIMPressureGradient( LevelSet2D ls, Grid2D gridP, Grid2D gridPX, Grid2D gridPY )
{
  PetscErrorCode ierr;
  IrregularNode *n;
  PetscReal s,
            **p  = gridP->v2,
            **px = gridPX->v2,
            **py = gridPY->v2,
            **phi= ls->g2d->v2;
  int i,j, d1=gridP->d1, d2=gridP->d2, x, y;
  
  PetscFunctionBegin;
  PetscLogEventBegin(EVENT_IIMPressureGradient,0,0,0,0);
  
  for( j = 0; j < d2; ++j) //px, py periodic bc
  {
    for( i = 0; i < d1; ++i)
    {
      px[i][j] = (p[(i+1)%d1][j] - p[(i-1+d1)%d1][j]) / 2.;
      py[i][j] = (p[i][(j+1)%d2] - p[i][(j-1+d2)%d2]) / 2.;
    }
  }
  
  for( i=0; i<ls->irregularNodes->len; i++)
  {
    n = &g_array_index( ls->irregularNodes, IrregularNode, i);
    x = n->x;
    y = n->y;
    s = n->sign;
    // X - gradient
    if( s == PetscSign( phi[x+1][y] ) )
    {
      if( s == PetscSign( phi[x-1][y] ) )
      {
        px[x][y] = (p[x+1][y] - p[x-1][y])/2.;//( )---( )---( )
      } else {
        px[x][y] = p[x+1][y] - p[x][y];      // (*)---( )---( )
      }
    } else {
      if( s == PetscSign( phi[x-1][y] ) )
      {
        px[x][y] = p[x][y] - p[x-1][y];      // ( )---( )---(*)
      } else {
        px[x][y] = IIMxGradient(n,p);        // (*)---( )---(*)
      }
    }
    // Y - gradient
    if( s == PetscSign( phi[x][y+1] ) )
    {
      if( s == PetscSign( phi[x][y-1] ) )
      {
        py[x][y] = (p[x][y+1] - p[x][y-1])/2.;
      } else {
        py[x][y] = p[x][y+1] - p[x][y];
      }
    } else {
      if( s == PetscSign( phi[x][y-1] ) )
      {
        py[x][y] = p[x][y] - p[x][y-1];
      } else {
        py[x][y] = IIMyGradient(n,p);
      }
    }
  }
  PetscLogEventEnd(EVENT_IIMPressureGradient,0,0,0,0);
  PetscFunctionReturn(0);
}

PetscReal IIMxGradient(IrregularNode *n, PetscReal **P )
{
  SETERRQ(PETSC_ERR_SUP, "IIMxGradient deprecated, use IIM Simple");
  PetscReal jp, jpx, jpy, p, pl, px,
    X = n->x + n->ox,     // Orthogonal projection  
    Y = n->y + n->oy,
    xi = n->x,            // Irregular Node
    yj = n->y,
    s = n->sign;
  PetscInt xl = xi + .5 < X  ? n->x - 1 : n->x + 1;
  
  p  = P[n->x][n->y];
  pl = P[xl][n->y];
  
  jp = n->f1;
  jpx = n->df2 * n->nx - n->df1 * n->ny;
  jpy = n->df2 * n->ny + n->df1 * n->nx;
  
  px  = p - pl - s*jp - s*jpx*(xl-X) - s*jpy*(yj-Y);
  px /= (xi - xl);
  return px * 0;
}

PetscReal IIMyGradient(IrregularNode *n, PetscReal **P )
{
  SETERRQ(PETSC_ERR_SUP, "IIMyGradient deprecated, use IIM Simple");
  PetscReal jp, jpx, jpy, p, pl, py,
    X = n->x + n->ox,     // Orthogonal projection  
    Y = n->y + n->oy,
    xi = n->x,            // Irregular Node
    yj = n->y,
    s = n->sign;
  PetscInt yl = yj + .5 < Y ? n->y - 1 : n->y + 1;
  
  jp = n->f1;
  jpx = n->df2 * n->nx - n->df1 * n->ny;
  jpy = n->df2 * n->ny + n->df1 * n->nx;
  
  p  = P[n->x][n->y];
  pl = P[n->x][yl];
  
  py  = p - pl - s*jp - s*jpx*(xi-X) - s*jpy*(yl-Y);
  py /= (yj - yl);
  return py * 0;
}

void rreg()
{
  PetscLogEventRegister(&EVENT_IIMPressureGradient,"IIMPressureGradient", 0);
}