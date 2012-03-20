#include "Grid.h"

PetscErrorCode Grid_MakeGrid( Grid g );
PetscErrorCode Grid_RestoreGrid( Grid g );

#undef __FUNCT__
#define __FUNCT__ "GridCreate"
PetscErrorCode GridCreate( Coor dh, iCoor pos, iCoor size, int dof, Grid *grid )
{
  Grid g;
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  ierr = PetscNew( struct _Grid, &g); CHKERRQ(ierr);
  g->d = dh;
  g->n = size;
  g->p = pos;
  g->dof = dof;
  ierr = PetscStrcpy(g->name,"grid"); CHKERRQ(ierr);
  
  g->is2D = size.z < 2 ? PETSC_TRUE : PETSC_FALSE; //TODO: alert user that for z < 2, will create a 2D grid
  g->SIZE = g->is2D ? size.x*size.y*dof : size.x*size.y*size.z*dof;
  g->MAXSIZE = g->SIZE;
  
  ierr = PetscMalloc(g->SIZE*sizeof(PetscReal), &g->v1); CHKERRQ(ierr);
  ierr = VecCreateSeqWithArray(PETSC_COMM_SELF, g->SIZE, g->v1, &g->v); CHKERRQ(ierr);
  ierr = VecZeroEntries(g->v); CHKERRQ(ierr); // TODO: is this redundant?
  ierr = Grid_MakeGrid(g); CHKERRQ(ierr);
  
  g->Interpolate = g->is2D ? GridInterpolate2D : GridInterpolate3D;

  *grid = g;
  
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "GridDestroy"
PetscErrorCode GridDestroy( Grid g )
{
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  if( g->filePos ) {
    ierr = PetscViewerDestroy(&g->filePos); CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&g->fileSize); CHKERRQ(ierr);
  }
  ierr = PetscFree(g->v1); CHKERRQ(ierr);
  ierr = VecDestroy(&g->v); CHKERRQ(ierr);
  ierr = PetscFree( g ); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "GridDuplicate"
PetscErrorCode GridDuplicate( Grid g, Grid *newg )
{
  Grid newgrid;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = GridCreate(g->d,g->p,g->n,g->dof,&newgrid); CHKERRQ(ierr);
  *newg = newgrid;
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "Grid_MakeGrid"
PetscErrorCode Grid_MakeGrid( Grid g )
{
  int dof = g->dof;
  iCoor s = g->n;
  iCoor p = g->p;
  void *grid=NULL;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  if( dof == 1 ) {
    if( g->is2D ) {
      ierr = VecGetArray2d( g->v, s.y, s.x, p.y, p.x,(PetscReal ***)&grid); CHKERRQ(ierr);
    } else {
      ierr = VecGetArray3d( g->v, s.z, s.y, s.x, p.z, p.y, p.x, (PetscReal ****)&grid); CHKERRQ(ierr);
    }
  } else {
    if( g->is2D ) {
      ierr = VecGetArray3d( g->v, s.y, s.x, dof, p.y, p.x, 0, (PetscReal ****)&grid); CHKERRQ(ierr);
    } else {
      ierr = VecGetArray4d( g->v, s.z, s.y, s.x, dof, p.z, p.y, p.x, 0, (PetscReal *****)&grid); CHKERRQ(ierr);
    }
  }
  g->grid = grid;
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "Grid_RestoreGrid"
PetscErrorCode Grid_RestoreGrid( Grid g )
{
  int dof = g->dof;
  iCoor s = g->n;
  iCoor p = g->p;
  void *grid = g->grid;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  if( dof == 1 ) {
    if( g->is2D ) {
      ierr = VecRestoreArray2d( g->v, s.y, s.x, p.y, p.x,(PetscReal ***)&grid); CHKERRQ(ierr);
    } else {
      ierr = VecRestoreArray3d( g->v, s.z, s.y, s.x, p.z, p.y, p.x, (PetscReal ****)&grid); CHKERRQ(ierr);
    }
  } else {
    if( g->is2D ) {
      ierr = VecRestoreArray3d( g->v, s.y, s.x, dof, p.y, p.x, 0, (PetscReal ****)&grid); CHKERRQ(ierr);
    } else {
      ierr = VecRestoreArray4d( g->v, s.z, s.y, s.x, dof, p.z, p.y, p.x, 0, (PetscReal *****)&grid); CHKERRQ(ierr);
    }
  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "GridResize"
PetscErrorCode GridResize( Grid g, iCoor pos, iCoor size )
{

  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = Grid_RestoreGrid(g); CHKERRQ(ierr);

  g->p = pos;
  g->n = size;
  g->SIZE = g->is2D ? size.x*size.y*g->dof : size.x*size.y*size.z*g->dof;

  //If requested size is smaller than previously allocated, enlarge
  if( g->SIZE > g->MAXSIZE )
  {
    ierr = PetscFree(g->v1); CHKERRQ(ierr);
    ierr = PetscMalloc(g->SIZE*sizeof(PetscReal), &g->v1); CHKERRQ(ierr);
    g->MAXSIZE = g->SIZE;
    ierr = PetscInfo3(0,"%s resizing: %d (%d MB)\n",g->name, g->MAXSIZE, g->MAXSIZE/(1024*1024) ); CHKERRQ(ierr);
  }
  ierr = VecDestroy(&g->v); CHKERRQ(ierr);
  ierr = VecCreateSeqWithArray(PETSC_COMM_SELF, g->SIZE, g->v1, &g->v); CHKERRQ(ierr);
  ierr = PetscMemzero(g->v1,g->MAXSIZE*sizeof(PetscReal)); CHKERRQ(ierr);
  ierr = Grid_MakeGrid(g); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode GridSetDx( Grid g, Coor d )
{
  g->d = d;
  return  0;
}

PetscErrorCode GridSetName( Grid g, const char *name )
{
  PetscErrorCode ierr;
  ierr = PetscStrcpy(g->name,name); CHKERRQ(ierr);
  if( g->filePos ) {
    ierr = PetscViewerDestroy(&g->filePos); CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&g->fileSize); CHKERRQ(ierr);
  }
  g->filePos = NULL;
  g->fileSize = NULL;
  return 0;
}

PetscErrorCode GridGet( Grid g, void *grid )
{
  *((void**)grid) = g->grid;
  return 0;
}

PetscErrorCode GridGetBounds( Grid g, iCoor *p, iCoor *q )
{
  *p = g->p;
  q->x = g->p.x + g->n.x;
  q->y = g->p.y + g->n.y;
  q->z = g->p.z + g->n.z;
  return 0;
}

double Bilinear2D( GridFunction2D gf, PetscReal **v2, Coor dh, double x, double y)
{
  int k,l;
  PetscReal xb, yb, sum = 0.;
  PetscInt xi = (int)floor(x);
  PetscInt yi = (int)floor(y);
  PetscReal xf = ( 2. * ( x - xi ) - 1. );
  PetscReal yf = ( 2. * ( y - yi ) - 1. );
  
  for( k = 0; k < 2; ++k)
  {
    for( l = 0; l < 2; ++l)
    {
      xb = 1. + (2. * k - 1.) * xf;
      yb = 1. + (2. * l - 1.) * yf;
      sum += gf( v2, xi+k, yi+l, dh ) * xb * yb;
    }
  }
  
  sum /= 4.;
  
  return sum;
}

inline double GridFunction2D_Identity( PetscReal **p, int i, int j, Coor d)
{
  return p[j][i];
}

inline double GridFunction2D_DerivX( PetscReal **p, int i, int j, Coor d)
{
  return (p[j][i+1]-p[j][i-1]) / (2. * d.x);
}

inline double GridFunction2D_DerivY( PetscReal **p, int i, int j, Coor d)
{
  return (p[j+1][i]-p[j-1][i]) / (2. * d.y);
}
/*Ex usage:
 *px = Bilinear2D( GridFunction2D_DerivX, ls->g2d, ox, oy ); 
 */
inline double GridFunction2D_Curv( double **p, int i, int j, Coor d )
{
  double px, py, px2, py2, pxx, pyy, pxy, k;
  px = (p[j][i+1]-p[j][i-1]) / (2. * d.x);
  py = (p[j+1][i]-p[j-1][i]) / (2. * d.y);
  px2 = px * px;
  py2 = py * py;
  pxx = (p[j][i-1] - 2.*p[j][i] + p[j][i+1]) / (d.x*d.x);
  pyy = (p[j-1][i] - 2.*p[j][i] + p[j+1][i]) / (d.y*d.y);
  pxy = (p[j-1][i-1] + p[j+1][i+1] - p[j+1][i-1] - p[j-1][i+1]) / (4.*d.x*d.y);
  k   = (pxx*py2 - 2.*px*py*pxy + pyy*px2) / sqrt((px2+py2)*(px2+py2)*(px2+py2));
/*
  const double TOL = 1e-1;
  const double mag = sqrt(px2+py2);
  if( PetscAbs(px/mag) < TOL || PetscAbs(py/mag) < TOL)
    k = 0;
*/
  if( k != k ) {
    //TODO: what to do when curvature is NaN? why?
    PetscInfo5(0,"Curvature NaN: px = %f, py = %f, pxx = %f, pyy = %f, pxy = %f\n",px,py,pxx,pyy,pxy);
    k = 0;
  }
  return k;
}

#undef __FUNCT__
#define __FUNCT__ "GridDrawBorder"
PetscErrorCode GridDrawBorder( Grid g, int borderWidth, PetscReal fill )
{
  int i,j,k;
  int b;
  iCoor p,q;
  PetscReal *phi=0;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = GridGet(g,&phi); CHKERRQ(ierr);
  if( g->is2D ) {
    PetscReal **phi2D = (PetscReal**)phi;
    for ( i = 0; i < g->n.x; ++i) {
      for (b = 0; b < borderWidth; ++b) {
        phi2D[         b][i] = fill;
        phi2D[g->n.y-1-b][i] = fill;
      }
    }
    for ( j = 0; j < g->n.y; ++j) {
      for (b = 0; b < borderWidth; ++b) {
        phi2D[j][         b] = fill;
        phi2D[j][g->n.x-1-b] = fill;
      }
    }
  } else {
    PetscReal ***phi3D = (PetscReal***)phi;
    phi3D[0][0][0] = 0;
    for ( k = 0; k < g->n.z; ++k) {
      for ( j = 0; j < g->n.y; ++j) {
        for ( b = 0; b < borderWidth; ++b) {
          phi3D[k][j][         b] = fill;
          phi3D[k][j][g->n.x-1-b] = fill;
        } // b
      } // j
      for ( b = 0; b < borderWidth; ++b) {
        for ( i = 0; i < g->n.x; ++i) {
          phi3D[k][         b][i] = fill;
          phi3D[k][g->n.y-1-b][i] = fill;
        } // i
      } // b
    } // k
    for ( b = 0; b < borderWidth; ++b) {
      for ( j = 0; j < g->n.y; ++j) {
        for ( i = 0; i < g->n.x; ++i) {
          phi3D[         b][j][i] = fill;
          phi3D[g->n.z-1-b][j][i] = fill;
        } // i
      } // j
    } // b
  } // 3D
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "GridFillRectangle"
PetscErrorCode GridFillRectangle( Grid g, Coor lo, Coor hi, PetscReal fill)
{
  Coor dh = g->d;
  iCoor p, q;
  int i,j,t;
  PetscReal **phi=0;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = GridGetBounds(g,&p,&q); CHKERRQ(ierr);
  t = lo.x/dh.x;
  p.x = t < p.x ? p.x : t;
  t = lo.y/dh.y;
  p.y = t < p.y ? p.y : t;

  t = hi.x/dh.x;
  q.x = q.x < t ? q.x : t;
  t = hi.y/dh.y;
  q.y = q.y < t ? q.y : t;

  ierr = GridGet(g,&phi); CHKERRQ(ierr);
  for (j = p.y; j < q.y; ++j) {
    for (i = p.x; i < q.x; ++i) {
      phi[j][i] = fill;
    }
  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "GridInterpolate"
PetscErrorCode GridInterpolate( Grid g, Coor X, PetscReal *val)
{
  PetscErrorCode ierr;
  PetscFunctionBegin;
  ierr = g->Interpolate(g,X,val); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "GridInterpolate2D"
PetscErrorCode GridInterpolate2D( Grid g, Coor X, PetscReal *val)
{
  int i,j;
  int xs, ys;
  PetscReal sx, sy;
  PetscReal **field=0;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = GridGet(g,&field); CHKERRQ(ierr);
  xs = (int)floor( X.x );
  ys = (int)floor( X.y );
  sx = X.x - xs;
  sy = X.y - ys;
  *val = 0;
  for( j = 0; j < 2; ++j)
  {
    for( i = 0; i < 2; ++i)
    {
      *val += ((1 - sx) * (1 - i) + i * sx) *
              ((1 - sy) * (1 - j) + j * sy) * field[ys+j][xs+i];
    }
  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "GridInterpolate3D"
PetscErrorCode GridInterpolate3D( Grid g, Coor X, PetscReal *val)
{
  int i,j,k;
  int xs, ys, zs;
  PetscReal sx, sy, sz;
  PetscReal ***field=0;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = GridGet(g,&field); CHKERRQ(ierr);
  xs = (int)floor( X.x );
  ys = (int)floor( X.y );
  zs = (int)floor( X.z );
  sx = X.x - xs;
  sy = X.y - ys;
  sz = X.z - zs;
  *val = 0;
  for( k = 0; k < 2; ++k)
  {
    for( j = 0; j < 2; ++j)
    {
      for( i = 0; i < 2; ++i)
      {
        *val += ((1 - sx) * (1 - i) + i * sx) *
                ((1 - sy) * (1 - j) + j * sy) *
                ((1 - sz) * (1 - k) + k * sz) * field[zs+k][ys+j][xs+i];
      }
    }
  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "GridCopy"
PetscErrorCode GridCopy( Grid g, Grid copy )
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = GridResize(copy,g->p,g->n); CHKERRQ(ierr);
  ierr = VecCopy(g->v,copy->v); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
