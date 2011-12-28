/*
 * Takes into account the staggered nature of the velocity grid
 *
 */
#include "Common.h"

PetscErrorCode InterpolateVelocity2D( const int udof, PetscReal ***field, const Coor X, Coor *vel )
{
  int dof,i,j;
  int xs, ys;
  PetscReal sx, sy;
  PetscReal sum;
  PetscReal *v = &(vel->x);

  for( dof = 0; dof < 2; dof++ )
  {
    xs = (int)floor(X.x + Tensor1[dof][0]);
    ys = (int)floor(X.y + Tensor1[dof][1]);
    sx = X.x - xs + Tensor1[dof][0];
    sy = X.y - ys + Tensor1[dof][1];
    sum = 0;
    for( j = 0; j < 2; ++j)
    {
      for( i = 0; i < 2; ++i)
      {
        sum += ((1 - sx) * (1 - i) + i * sx) *
               ((1 - sy) * (1 - j) + j * sy) * field[ys+j][xs+i][dof+udof];
      }
    }
    v[dof] = sum;
  }
  return 0;
}

PetscErrorCode InterpolateVelocity3D(const int udof, PetscReal ****field, const Coor X, Coor *vel )
{
  int dof, i,j,k;
  int xs,ys,zs;
  PetscReal sx,sy,sz;
  PetscReal sum;
  PetscReal *v = &(vel->x);

  for( dof = 0; dof < 3; dof++ )
  {
    xs = (int)floor(X.x + Tensor1[dof][0]);
    ys = (int)floor(X.y + Tensor1[dof][1]);
    zs = (int)floor(X.z + Tensor1[dof][2]);
    sx = X.x - xs + Tensor1[dof][0];
    sy = X.y - ys + Tensor1[dof][1];
    sz = X.z - zs + Tensor1[dof][2];
    sum = 0;
    for( k = 0; k < 2; ++k)
    {
      for( j = 0; j < 2; ++j)
      {
        for( i = 0; i < 2; ++i)
        {
          sum += ((1 - sx) * (1 - i) + i * sx) *
                 ((1 - sy) * (1 - j) + j * sy) *
                 ((1 - sz) * (1 - k) + k * sz) * field[zs+k][ys+j][xs+i][dof+udof];
        }
      }
    }
    v[dof] = sum;
  }
  return 0;
}
