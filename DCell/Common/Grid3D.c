#include "Grid.h"

// Equ 4.36 - 4.37 
double Bilinear3D( GridFunction3D gf, PetscReal ***v3, Coor dh, Coor p )
{
  int i,j,k;
  PetscReal xb, yb, zb, sum = 0.0;
  PetscInt xs = floor(p.x);
  PetscInt ys = floor(p.y);
  PetscInt zs = floor(p.z);
  PetscReal sx = p.x - xs,
            sy = p.y - ys,
            sz = p.z - zs;

  
  for( i = 0; i < 2; ++i)
  {
    for( j = 0; j < 2; ++j)
    {
      for( k = 0; k < 2; ++k)
      {
        xb = i * sx + (1-i) * (1-sx);
        yb = j * sy + (1-j) * (1-sy);
        zb = k * sz + (1-k) * (1-sz);
        sum += gf( v3, xs+i, ys+j, zs+k, dh ) * xb * yb * zb;
      }
    }
  }
  
  sum /= 8.;
  
  return sum;
}

inline double GridFunction3D_Identity( PetscReal ***p, int i, int j, int k, Coor d)
{
  return p[k][j][i];
}

inline double GridFunction3D_DerivX( PetscReal ***p, int i, int j, int k, Coor d)
{
  return (p[k][j][i+1]-p[k][j][i-1]) / (2.0 * d.x);
}

inline double GridFunction3D_DerivY( PetscReal ***p, int i, int j, int k, Coor d)
{
  return (p[k][j+1][i]-p[k][j-1][i]) / (2.0 * d.y);
}

inline double GridFunction3D_DerivZ( PetscReal ***p, int i, int j, int k, Coor d)
{
  return (p[k+1][j][i]-p[k-1][j][i]) / (2.0 * d.z);
}

inline double GridFunction3D_DerivXX( PetscReal ***p, int i, int j, int k, Coor d)
{
  return (p[k][j][i-1] - 2.*p[k][j][i] + p[k][j][i+1]) / (d.x*d.x);
}

inline double GridFunction3D_DerivYY( PetscReal ***p, int i, int j, int k, Coor d)
{
  return (p[k][j-1][i] - 2.*p[k][j][i] + p[k][j+1][i]) / (d.y*d.y);
}

inline double GridFunction3D_DerivZZ( PetscReal ***p, int i, int j, int k, Coor d)
{
  return (p[k-1][j][i] - 2.*p[k][j][i] + p[k+1][j][i]) / (d.z*d.z);
}

inline double GridFunction3D_DerivXY( PetscReal ***p, int i, int j, int k, Coor d)
{
  return (p[k][j-1][i-1] + p[k][j+1][i+1] - p[k][j+1][i-1] - p[k][j-1][i+1] ) / (4.*d.x*d.y);
}
  
inline double GridFunction3D_DerivXZ( PetscReal ***p, int i, int j, int k, Coor d)
{
  return (p[k-1][j][i-1] + p[k+1][j][i+1] - p[k+1][j][i-1] - p[k-1][j][i+1] ) / (4.*d.x*d.z);
}

inline double GridFunction3D_DerivYZ( PetscReal ***p, int i, int j, int k, Coor d)
{
  return (p[k-1][j-1][i] + p[k+1][j+1][i] - p[k+1][j-1][i] - p[k-1][j+1][i] ) / (4.*d.y*d.z);
}

inline double GridFunction3D_Curv( double ***p, int i, int j, int k, Coor d )
{
  double px, py, pz, px2, py2, pz2, pxx, pyy, pzz,
         ppp, pxy, pxz, pyz, c;
  px = GridFunction3D_DerivX( p, i, j, k, d);
  py = GridFunction3D_DerivY( p, i, j, k, d);
  pz = GridFunction3D_DerivZ( p, i, j, k, d);
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
  return c;
}
