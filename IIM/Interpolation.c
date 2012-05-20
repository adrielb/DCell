#include "ImmersedInterfaceMethod.h"

inline PetscReal IIMVelocityCorrection(
    PetscReal x, PetscReal ux1, PetscReal a1, PetscReal ux2, PetscReal a2,
    PetscReal y, PetscReal uy1, PetscReal uy2 )
{
  const PetscReal C12 = x < a1 ? (1-a1)*x*ux1 : a1*(1-x)*ux1;
  const PetscReal C34 = x < a2 ? (1-a2)*x*ux2 : a2*(1-x)*ux2;
  const PetscReal C13 = y < b1 ? (1-b1)*y*uy1 : b1*(1-y)*uy1;
  const PetscReal C24 = y < b2 ? (1-b2)*y*uy2 : b2*(1-y)*uy2;

  return (1-y) * C12 + y * C34 +
         (1-x) * C13 + x * C24;
}

inline PetscReal IIMVelocityCorrection( PetscReal x, PetscReal y, PetscReal z, PetscReal ux1, PetscReal a1 )
{
  const PetscReal C12 = X.x < a1 ? (1-a1) * X.x * ux1 : a1 * (1-X.x) * ux1;

  return (1-y) * C12;
}

PetscErrorCode InterpolateVelocity2D( const int udof, const Coor X, Coor *vel )
{
  int dof,i,j;
  int xs, ys;
  PetscReal sx, sy;
  PetscReal sum;
  PetscReal *v = &(vel->x);

  IrregularNode *n;

  n->

  for( dof = 0; dof < 2; dof++ )
  {
    xs = (int)floor(X.x + Tensor1[dof][0]);
    ys = (int)floor(X.y + Tensor1[dof][1]);
    sx = X.x - xs + Tensor1[dof][0];
    sy = X.y - ys + Tensor1[dof][1];
    sum = 0;
    for( j = 0; j < 2; ++j) {
      for( i = 0; i < 2; ++i) {
        sum += ((1 - sx) * (1 - i) + i * sx) *
               ((1 - sy) * (1 - j) + j * sy) * field[ys+j][xs+i][dof+udof];
      }
    }
    v[dof] = sum;
    v[dof] += IIMVelocityCorrection(
        sx, ujx[ys][xs], ujx[ys+1][xs],
        sy, ujy[ys][xs], ujy[ys][xs+1] );
  }
  return 0;
}
