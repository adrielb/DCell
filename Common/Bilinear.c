#include "Grid.h"

double GridBilinear( Grid g, GridFunction gf, Coor p )
{
  int i,j,k;
  PetscReal *val;
  PetscReal xb, yb, zb;
  PetscReal sum = 0.0;
  iCoor S;
  Coor s;
  CoorToIndex2(g->aabb.lo, g->d, p, &S, &s);
  GridGet( g, &val );

  s.x -= S.x;
  s.y -= S.y;
  s.z -= S.z;

  for( i = 0; i < 2; ++i)
  {
    for( j = 0; j < 2; ++j)
    {
      for( k = 0; k < 2; ++k)
      {
        xb = i * s.x + (1-i) * (1-s.x);
        yb = j * s.y + (1-j) * (1-s.y);
        zb = k * s.z + (1-k) * (1-s.z);
        sum += gf( val, (iCoor){S.x+i, S.y+j, S.z+k}, g->d ) * xb * yb * zb;
      }
    }
  }
  
  return sum;
}

inline double GridFunction2D_DerivX( PetscReal *v, iCoor p, Coor d)
{
  PetscReal **v2 = (PetscReal**)v;
  return (v2[p.y][p.x+1] - v2[p.y][p.x-1]) / (2. * d.x);
}

inline double GridFunction2D_DerivY( PetscReal *v, iCoor p, Coor d)
{
  PetscReal **v2 = (PetscReal**)v;
  return (v2[p.y+1][p.x] - v2[p.y-1][p.x]) / (2. * d.y);
}

inline double GridFunction2D_Curv( PetscReal *v, iCoor x, Coor d )
{
  PetscReal **p = (PetscReal**)v;
  const int i = x.x;
  const int j = x.y;
  double px, py, px2, py2, pxx, pyy, pxy, k;
  px = (p[j][i+1]-p[j][i-1]) / (2. * d.x);
  py = (p[j+1][i]-p[j-1][i]) / (2. * d.y);
  px2 = px * px;
  py2 = py * py;
  pxx = (p[j][i-1] - 2.*p[j][i] + p[j][i+1]) / (d.x*d.x);
  pyy = (p[j-1][i] - 2.*p[j][i] + p[j+1][i]) / (d.y*d.y);
  pxy = (p[j-1][i-1] + p[j+1][i+1] - p[j+1][i-1] - p[j-1][i+1]) / (4.*d.x*d.y);
  k   = (pxx*py2 - 2.*px*py*pxy + pyy*px2) / sqrt((px2+py2)*(px2+py2)*(px2+py2));

  if( k != k ) {
    //TODO: what to do when curvature is NaN? why?
    PetscInfo5(0,"Curvature NaN: px = %f, py = %f, pxx = %f, pyy = %f, pxy = %f\n",px,py,pxx,pyy,pxy);
    k = 0;
  }
  return k;
}

inline double GridFunction3D_DerivX( PetscReal *v, iCoor x, Coor d)
{
  PetscReal ***p = (PetscReal***)v;
  return (p[x.z][x.y][x.x+1]-p[x.z][x.y][x.x-1]) / (2.0 * d.x);
}

inline double GridFunction3D_DerivY( PetscReal *v, iCoor x, Coor d)
{
  PetscReal ***p = (PetscReal***)v;
  return (p[x.z][x.y+1][x.x]-p[x.z][x.y-1][x.x]) / (2.0 * d.y);
}

inline double GridFunction3D_DerivZ( PetscReal *v, iCoor x, Coor d)
{
  PetscReal ***p = (PetscReal***)v;
  return (p[x.z+1][x.y][x.x]-p[x.z-1][x.y][x.x]) / (2.0 * d.z);
}

inline double GridFunction3D_DerivXX( PetscReal *v, iCoor x, Coor d)
{
  PetscReal ***p = (PetscReal***)v;
  return (p[x.z][x.y][x.x-1] - 2.*p[x.z][x.y][x.x] + p[x.z][x.y][x.x+1]) / (d.x*d.x);
}

inline double GridFunction3D_DerivYY( PetscReal *v, iCoor x, Coor d)
{
  PetscReal ***p = (PetscReal***)v;
  return (p[x.z][x.y-1][x.x] - 2.*p[x.z][x.y][x.x] + p[x.z][x.y+1][x.x]) / (d.y*d.y);
}

inline double GridFunction3D_DerivZZ( PetscReal *v, iCoor x, Coor d)
{
  PetscReal ***p = (PetscReal***)v;
  return (p[x.z-1][x.y][x.x] - 2.*p[x.z][x.y][x.x] + p[x.z+1][x.y][x.x]) / (d.z*d.z);
}

inline double GridFunction3D_DerivXY( PetscReal *v, iCoor x, Coor d)
{
  PetscReal ***p = (PetscReal***)v;
  return (p[x.z][x.y-1][x.x-1] + p[x.z][x.y+1][x.x+1] - p[x.z][x.y+1][x.x-1] - p[x.z][x.y-1][x.x+1] ) / (4.*d.x*d.y);
}
  
inline double GridFunction3D_DerivXZ( PetscReal *v, iCoor x, Coor d)
{
  PetscReal ***p = (PetscReal***)v;
  return (p[x.z-1][x.y][x.x-1] + p[x.z+1][x.y][x.x+1] - p[x.z+1][x.y][x.x-1] - p[x.z-1][x.y][x.x+1] ) / (4.*d.x*d.z);
}

inline double GridFunction3D_DerivYZ( PetscReal *v, iCoor x, Coor d)
{
  PetscReal ***p = (PetscReal***)v;
  return (p[x.z-1][x.y-1][x.x] + p[x.z+1][x.y+1][x.x] - p[x.z+1][x.y-1][x.x] - p[x.z-1][x.y+1][x.x] ) / (4.*d.y*d.z);
}

inline double GridFunction3D_Curv( PetscReal *v, iCoor x, Coor d)
{
  PetscReal ***p = (PetscReal***)v;
  const int i = x.x;
  const int j = x.y;
  const int k = x.z;
  double px, py, pz, px2, py2, pz2, pxx, pyy, pzz,
         ppp, pxy, pxz, pyz, c;
  px = GridFunction3D_DerivX( v, x, d);
  py = GridFunction3D_DerivY( v, x, d);
  pz = GridFunction3D_DerivZ( v, x, d);
  px2 = px * px;
  py2 = py * py;
  pz2 = pz * pz;
  ppp = px2+py2+pz2;
  pxx = (p[k][j][i-1] - 2.*p[k][j][i] + p[k][j][i+1]) / (d.x*d.x);
  pyy = (p[k][j-1][i] - 2.*p[k][j][i] + p[k][j+1][i]) / (d.y*d.y);
  pzz = (p[k-1][j][i] - 2.*p[k][j][i] + p[k+1][j][i]) / (d.z*d.z);
  pxy = (p[k][j-1][i-1] + p[k][j+1][i+1] - p[k][j+1][i-1] - p[k][j-1][i+1] ) / (4.*d.x*d.y);
  pxz = (p[k-1][j][i-1] + p[k+1][j][i+1] - p[k+1][j][i-1] - p[k-1][j][i+1] ) / (4.*d.x*d.z);
  pyz = (p[k-1][j-1][i] + p[k+1][j+1][i] - p[k+1][j-1][i] - p[k-1][j+1][i] ) / (4.*d.y*d.z);
  c   = px2*(pyy+pzz) + py2*(pxx+pzz) + pz2*(pxx+pyy);
  c  -= 2.*(px*py*pxy + px*pz*pxz + py*pz*pyz);
  c  /= sqrt(ppp*ppp*ppp);

  if( c != c ) {
    c = 0;
  }
  return c;
}

#undef __FUNCT__
#define __FUNCT__ "GridFillBox"
PetscErrorCode GridFillBox( Grid g, Coor lo, Coor hi, PetscReal fill)
{
  Coor dh = g->d;
  iCoor p, q;
  int i,j,k,t;
  PetscReal *phi=0;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = GridGetBounds(g,&p,&q); CHKERRQ(ierr);
  t = (int)(lo.x/dh.x);
  p.x = t < p.x ? p.x : t;
  t = (int)(lo.y/dh.y);
  p.y = t < p.y ? p.y : t;
  t = (int)(lo.z/dh.z);
  p.z = t < p.z ? p.z : t;

  t = (int)(hi.x/dh.x);
  q.x = q.x < t ? q.x : t;
  t = (int)(hi.y/dh.y);
  q.y = q.y < t ? q.y : t;
  t = (int)(hi.z/dh.z);
  q.z = q.z < t ? q.z : t;

  ierr = GridGet(g,&phi); CHKERRQ(ierr);
  if( g->is2D ) {
    PetscReal **phi2D = (PetscReal**)phi;
    for (j = p.y; j < q.y; ++j) {
      for (i = p.x; i < q.x; ++i) {
        phi2D[j][i] = fill;
      }
    }
  } else {
    PetscReal ***phi3D = (PetscReal***)phi;
    for (k = p.z; k < q.z; ++k) {
      for (j = p.y; j < q.y; ++j) {
        for (i = p.x; i < q.x; ++i) {
          phi3D[k][j][i] = fill;
        }
      }
    }
  }
  PetscFunctionReturn(0);
}
