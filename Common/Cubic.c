#include "Grid.h"

#undef __FUNCT__
#define __FUNCT__ "GridInterpolate_Cubic"
PetscErrorCode GridInterpolate_Cubic( Grid g, Coor X, PetscReal *val)
{
  int i, i1, i2, i3;
  int j, j1, j2, j3;
  iCoor p,q;
  Coor s;   //
  iCoor S;  // lo of cell containing X
  PetscReal  **phi=0;
  PetscReal a[16];
  PetscReal x, x2, x3;
  PetscReal y, y2, y3;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = GridGet(g,&phi); CHKERRQ(ierr);
  GridGetBounds( g, &p, &q);
  q.x--;
  q.y--;
  q.z--;
  CoorToIndex2(g->aabb.lo, g->d, X, &S, &s);

  i = S.x - 1;
  j = S.y - 1;
  x = s.x - S.x;
  y = s.y - S.y;
  x2 = x*x;
  x3 = x2*x;
  y2 = y*y;
  y3 = y2*y;

  i1 = i+1; if( i1 < p.x ) i1 = p.x;  if( i1 > q.x ) i1 = q.x;
  i2 = i+2; if( i2 < p.x ) i2 = p.x;  if( i2 > q.x ) i2 = q.x;
  i3 = i+3; if( i3 < p.x ) i3 = p.x;  if( i3 > q.x ) i3 = q.x;
            if( i  < p.x ) i  = p.x;  if( i  > q.x ) i  = q.x;


  j1 = j+1; if( j1 < p.y ) j1 = p.y;  if( j1 > q.y ) j1 = q.y;
  j2 = j+2; if( j2 < p.y ) j2 = p.y;  if( j2 > q.y ) j2 = q.y;
  j3 = j+3; if( j3 < p.y ) j3 = p.y;  if( j3 > q.y ) j3 = q.y;
            if( j  < p.y ) j  = p.y;  if( j  > q.y ) j  = q.y;

  a[0] = phi[j1][i1];
  a[1] = (-phi[j1][i] + phi[j1][i2])/2.;
  a[2] = phi[j1][i] - 3*phi[j1][i1] + 2*phi[j1][i2] + (phi[j1][i1] - phi[j1][i3])/2.;
  a[3] = 2*phi[j1][i1] - 2*phi[j1][i2] + (-phi[j1][i] + phi[j1][i2])/2. + (-phi[j1][i1] + phi[j1][i3])/2.;
  a[4] = (-phi[j][i1] + phi[j2][i1])/2.;
  a[5] = (phi[j][i] - phi[j][i2] - phi[j2][i] + phi[j2][i2])/4.;
  a[6] = (-3*(-phi[j][i1] + phi[j2][i1]))/2. + (-phi[j][i] + phi[j][i2] + phi[j2][i] - phi[j2][i2])/2. + (3*(-phi[j][i2] + phi[j2][i2]))/2. + (-phi[j][i1] + phi[j][i3] + phi[j2][i1] - phi[j2][i3])/4.;
  a[7] = -phi[j][i1] + phi[j][i2] + phi[j2][i1] - phi[j2][i2] + (phi[j][i] - phi[j][i2] - phi[j2][i] + phi[j2][i2])/4. + (phi[j][i1] - phi[j][i3] - phi[j2][i1] + phi[j2][i3])/4.;
  a[8] = phi[j][i1] - 3*phi[j1][i1] + 2*phi[j2][i1] + (phi[j1][i1] - phi[j3][i1])/2.;
  a[9] = (-3*(-phi[j1][i] + phi[j1][i2]))/2. + (-phi[j][i] + phi[j][i2] + phi[j2][i] - phi[j2][i2])/2. + (3*(-phi[j2][i] + phi[j2][i2]))/2. + (-phi[j1][i] + phi[j1][i2] + phi[j3][i] - phi[j3][i2])/4.;
  a[10] = phi[j][i] - phi[j][i2] + 9*phi[j1][i1] - 9*phi[j1][i2] + 3*(-phi[j1][i] + phi[j1][i2]) + (3*(-phi[j1][i1] + phi[j1][i3]))/2. - phi[j2][i] - 9*phi[j2][i1] + 3*(-phi[j][i1] + phi[j2][i1]) + 10*phi[j2][i2] - 3*(-phi[j][i2] + phi[j2][i2]) - 3*(-phi[j2][i] + phi[j2][i2]) - (3*(-phi[j2][i1] + phi[j2][i3]))/2. + (phi[j][i1] - phi[j][i3] - phi[j2][i1] + phi[j2][i3])/2. + (3*(-phi[j1][i1] + phi[j3][i1]))/2. - (3*(-phi[j1][i2] + phi[j3][i2]))/2. + (phi[j1][i] - phi[j1][i2] - phi[j3][i] + phi[j3][i2])/2. + (phi[j1][i1] - phi[j1][i3] - phi[j3][i1] + phi[j3][i3])/4.;
  a[11] = -5*phi[j1][i1] + 5*phi[j1][i2] - (3*(-phi[j1][i] + phi[j1][i2]))/2. - (3*(-phi[j1][i1] + phi[j1][i3]))/2. + 6*phi[j2][i1] - 2*(-phi[j][i1] + phi[j2][i1]) + (-phi[j][i] + phi[j][i2] + phi[j2][i] - phi[j2][i2])/2. - 6*phi[j2][i2] + 2*(-phi[j][i2] + phi[j2][i2]) + (3*(-phi[j2][i] + phi[j2][i2]))/2. + (-phi[j][i1] + phi[j][i3] + phi[j2][i1] - phi[j2][i3])/2. + (3*(-phi[j2][i1] + phi[j2][i3]))/2. - phi[j3][i1] + (-phi[j1][i] + phi[j1][i2] + phi[j3][i] - phi[j3][i2])/4. + phi[j3][i2] + (-phi[j1][i1] + phi[j1][i3] + phi[j3][i1] - phi[j3][i3])/4.;
  a[12] = 2*phi[j1][i1] - 2*phi[j2][i1] + (-phi[j][i1] + phi[j2][i1])/2. + (-phi[j1][i1] + phi[j3][i1])/2.;
  a[13] = -phi[j1][i] + phi[j1][i2] + phi[j2][i] - phi[j2][i2] + (phi[j][i] - phi[j][i2] - phi[j2][i] + phi[j2][i2])/4. + (phi[j1][i] - phi[j1][i2] - phi[j3][i] + phi[j3][i2])/4.;
  a[14] = -5*phi[j1][i1] + 6*phi[j1][i2] - 2*(-phi[j1][i] + phi[j1][i2]) - phi[j1][i3] + 5*phi[j2][i1] - (3*(-phi[j][i1] + phi[j2][i1]))/2. + (-phi[j][i] + phi[j][i2] + phi[j2][i] - phi[j2][i2])/2. - 6*phi[j2][i2] + (3*(-phi[j][i2] + phi[j2][i2]))/2. + 2*(-phi[j2][i] + phi[j2][i2]) + (-phi[j][i1] + phi[j][i3] + phi[j2][i1] - phi[j2][i3])/4. + phi[j2][i3] - (3*(-phi[j1][i1] + phi[j3][i1]))/2. + (-phi[j1][i] + phi[j1][i2] + phi[j3][i] - phi[j3][i2])/2. + (3*(-phi[j1][i2] + phi[j3][i2]))/2. + (-phi[j1][i1] + phi[j1][i3] + phi[j3][i1] - phi[j3][i3])/4.;
  a[15] = -phi[j][i1] + phi[j][i2] - phi[j1][i] + 2*phi[j1][i1] - 2*phi[j1][i2] + phi[j1][i3] + phi[j2][i] - 2*phi[j2][i1] + 2*phi[j2][i2] + (phi[j][i] - phi[j][i2] - phi[j2][i] + phi[j2][i2])/4. - phi[j2][i3] + (phi[j][i1] - phi[j][i3] - phi[j2][i1] + phi[j2][i3])/4. + phi[j3][i1] - phi[j3][i2] + (phi[j1][i] - phi[j1][i2] - phi[j3][i] + phi[j3][i2])/4. + (phi[j1][i1] - phi[j1][i3] - phi[j3][i1] + phi[j3][i3])/4.;

  *val = a[0] + x*a[1] + x2*a[2] + x3*a[3] + y*a[4] + x*y*a[5] + x2*y*a[6] + x3*y*a[7] + y2*a[8] + x*y2*a[9] + x2*y2*a[10] + x3*y2*a[11] + y3*a[12] + x*y3*a[13] + x2*y3*a[14] + x3*y3*a[15];

  PetscFunctionReturn(0);
}
