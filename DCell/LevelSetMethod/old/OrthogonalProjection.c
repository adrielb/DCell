#include "LevelSetMethod.h"

/* Orthogonal Projection solves:
 * d = Dp(xi) / | Dp(xi) |
 * p(xi) + | Dp(xi) | a + 0.5 dT.He.d a^2 == 0
 * X* = xi + a d
 * p(X*) == 0
 *
 * TODO: Need a more robust root finding procedure where the second order method switches to first order
 *       when the Hessian is singular
 */

#undef __FUNCT__
#define __FUNCT__ "OrthogonalProjection2D_1st"
PetscErrorCode OrthogonalProjection2D_1st( double phi[3][3], double PHI_INF, Coor *op )
{
  int m, i, j, x, y;
  double p, px, py, pxy, pn, a, f, df;
  double p0, nx,ny, hxx, hxy, hyy, pHp;
  Coor Xs;
  double mindist = 2, dist;
  const int nei[4][2] = {{-1,0},{0,-1},{1,0},{0,1}};
  
  PetscFunctionBegin;
  for (i = 0; i < 4; ++i) {
    x = nei[i][0] + 1;
    y = nei[i][1] + 1;
    if( phi[1][1] * phi[y][x] <= 0 ) {
      dist = phi[1][1] / (phi[1][1] - phi[y][x]);
      if( dist < mindist ) {
        mindist = dist;
        op->x = dist * nei[i][0];
        op->y = dist * nei[i][1];
      }
    }
  }
  PetscFunctionReturn(0);
  for (j = 0; j < 2; ++j) {
    for (i = 0; i < 2; ++i) {
      p   = phi[j][i];
      px  = phi[j][i+1] - phi[j][i];
      py  = phi[j+1][i] - phi[j][i];
      pxy = phi[j][i] + phi[j+1][i+1] - phi[j+1][i] - phi[j][i+1];
      x = 1 - i;
      y = 1 - j;
//      for (y = 0; y < 2; y++) {
//        for (x = 0; x < 2; x++) {
          p0 = p + px*x + py*y + pxy*x*y;
          nx = px + pxy*y;
          ny = py + pxy*x;
          pn = sqrt( nx*nx + ny*ny );
          hxx = 0;
          hxy = pxy;
          hyy = 0;
          pHp = hxx*nx*nx + 2*hxy*nx*ny + hyy*ny*ny / (nx*nx+ny*ny);

          // Newton's iteration
          a = 0;
          for ( m = 0; m < 5; ++m) {
            f = p0 + pn * a + 0.5*a*a*pHp;
            df = pn + a*pHp;
            a = a - f / df;
          }
          if( a != a ) continue;

          Xs.x = x + a * px / pn;
          Xs.y = y + a * py / pn;
          if( Xs.x < 0 || Xs.x > 1 || Xs.y < 0 || Xs.y > 1  ) continue;

          Xs.x = Xs.x + i - 1;
          Xs.y = Xs.y + j - 1;
          dist = sqrt(Xs.x*Xs.x + Xs.y*Xs.y);

          if( dist < mindist ) {
            mindist = dist;
            op->x = Xs.x;
            op->y = Xs.y;
          }
//        } // x
//      } // y
    } // i
  } // j

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "OrthogonalProjection2D"
PetscErrorCode OrthogonalProjection2D( double phi[3][3], double PHI_INF, Coor *op)
{
  int m;
  const int i = 1, j = 1; //makes indexing easier
  double p, px, py, px2, py2, pxx, pyy, pxy, pn, c, a, f, df;

  PetscFunctionBegin;

  p  = phi[j][i];
  px = ( phi[j][i+1] - phi[j][i-1] ) / 2.;
  py = ( phi[j+1][i] - phi[j-1][i] ) / 2.;
  px2 = px * px;
  py2 = py * py;
  pxx= ( phi[j][i+1] - 2.*p + phi[j][i-1] );
  pyy= ( phi[j+1][i] - 2.*p + phi[j-1][i] );
  pxy= ( phi[j-1][i-1] + phi[j+1][i+1] - phi[j+1][i-1] - phi[j-1][i+1] ) / 4.;
  pn = sqrt( px2 + py2 );

  c = ( px2 * pxx + 2 * px * py * pxy + py2 * pyy ) / ( px2+py2 );

  // Newton's iteration
  a = 0;
  for ( m = 0; m < 5; ++m) {
    f = p + pn * a + 0.5*a*a*c;
    df = pn + a*c;
    a = a - f / df;
  }

  if( PetscAbs(a) < 1 ) {
    op->x = a * px / pn;
    op->y = a * py / pn;
  } else {
    OrthogonalProjection2D_1st( phi, PHI_INF, op );
  }

  PetscFunctionReturn(0);
}

/*
#undef __FUNCT__
#define __FUNCT__ "OrthogonalProjection2D_2nd"
PetscErrorCode OrthogonalProjection2D_2nd( double phi[3][3], double PHI_INF, double *ox, double *oy  )
{
  int m;
  const int i = 1, j = 1; //makes indexing easier
  double p, px, py, pxx, pyy, pxy, pxxy, pxyy, pxxyy, pn, a, f, df;
  double x,y, p0, nx,ny, hxx, hxy, hyy, pHp;
  Coor Xs;
  double mindist = 2, dist;

  PetscFunctionBegin;

  p   = phi[j][i];
  px  = ( phi[j][i+1] - phi[j][i-1] ) / 2.;
  py  = ( phi[j+1][i] - phi[j-1][i] ) / 2.;
  pxx = ( phi[j][i+1] - 2.*p + phi[j][i-1] ) / 2.;
  pyy = ( phi[j+1][i] - 2.*p + phi[j-1][i] ) / 2.;
  pxy = ( phi[j-1][i-1] + phi[j+1][i+1] - phi[j+1][i-1] - phi[j-1][i+1] ) / 4.;
  pxxy = (-phi[-1 + j][-1 + i] + 2*phi[-1 + j][i] - phi[-1 + j][1 + i] + phi[1 + j][-1 + i] - 2*phi[1 + j][i] + phi[1 + j][1 + i]) / 4;
  pxyy = (-phi[-1 + j][-1 + i] + phi[-1 + j][1 + i] + 2*phi[j][-1 + i] - 2*phi[j][1 + i] - phi[1 + j][-1 + i] + phi[1 + j][1 + i]) / 4;
  pxxyy = (phi[-1 + j][-1 + i] - 2*phi[-1 + j][i] + phi[-1 + j][1 + i] - 2*phi[j][-1 + i] + 4*phi[j][i] - 2*phi[j][1 + i] + phi[1 + j][-1 + i] - 2*phi[1 + j][i] + phi[1 + j][1 + i]) / 4;

  for (y = -1; y < 1.1; y+=0.5) {
    for (x = -1; x < 1.1; x+=0.5) {
      p0 = p + px*x + py*y + pxy*x*y + pxx*x*x + pyy*y*y + pxxy*x*x*y + pxyy*x*y*y + pxxyy*x*x*y*y;
      nx = px + 2*pxx*x + pxy*y + 2*pxxy*x*y + pxyy*y*y + 2*pxxyy*x*y*y;
      ny = py + pxy*x + pxxy*x*x + 2*pyy*y + 2*pxyy*x*y + 2*pxxyy*x*x*y;
      pn = sqrt( nx*nx + ny*ny );
      hxx = 2*pxx + 2*pxxy*y + 2*pxxyy*y*y;
      hxy = pxy + 2*pxxy*x + 2*pxyy*y + 4*pxxyy*x*y;
      hyy = 2*pyy + 2*pxyy*x + 2*pxxyy*x*x;
      pHp = hxx*nx*nx + 2*hxy*nx*ny + hyy*ny*ny / (nx*nx+ny*ny);

      // Newton's iteration
      a = 0;
      for ( m = 0; m < 5; ++m) {
        f = p0 + pn * a + 0.5*a*a*pHp;
        df = pn + a*pHp;
        a = a - f / df;
      }
      if( a != a ) continue;

      Xs.x = x + a * px / pn;
      Xs.y = y + a * py / pn;
      dist = sqrt(Xs.x*Xs.x + Xs.y*Xs.y);

      if( dist < mindist ) {
        mindist = dist;
        *ox = Xs.x;
        *oy = Xs.y;
      }
    } // x
  } // y

  PetscFunctionReturn(0);
}



#undef __FUNCT__
#define __FUNCT__ "OrthogonalProjection2D"
PetscErrorCode OrthogonalProjection2D( double phi[3][3], double PHI_INF, double *ox, double *oy  )
{
  PetscErrorCode ierr;
  const int i = 1, j = 1; //makes indexing easier
  double p, px, py, px2, py2, pxx, pyy, pxy, pn, c, b, discr, a1, a2, a;

  PetscFunctionBegin;

  p  = phi[j][i];
  px = ( phi[j][i+1] - phi[j][i-1] ) / 2.;
  py = ( phi[j+1][i] - phi[j-1][i] ) / 2.;
  px2 = px * px;
  py2 = py * py;
  pxx= ( phi[j][i+1] - 2.*p + phi[j][i-1] );
  pyy= ( phi[j+1][i] - 2.*p + phi[j-1][i] );
  pxy= ( phi[j-1][i-1] + phi[j+1][i+1] - phi[j+1][i-1] - phi[j-1][i+1] ) / 4.;
  pn = sqrt( px2 + py2 );

  b = pn;
  c = ( px2 * pxx + 2 * px * py * pxy + py2 * pyy ) / ( 2 * ( px2+py2 ) );
  discr = b*b - 4 * p * c;

  for (int m = 0; m < 5; ++m) {
    f = p + pn * a + 0.5 a*a;
    df = pn + a;
    a = a - f / df;
  }


  if( c == 0.0 ) {
    a = -p / b;
    *ox = a * px / pn;
    *oy = a * py / pn;
    PetscFunctionReturn(0);
  }

  if( discr < 0 ) {
    ierr = OrthoProjFirstOrder2D(phi,PHI_INF,ox,oy); CHKERRQ(ierr);
    printf("discr: %e  c: %e {%f,%f}\n", discr, c, *ox, *oy);

  }

  discr = sqrt(discr);
  a1 = ( -b + discr ) / (2*c);
  a2 = ( -b - discr ) / (2*c);
  a = PetscAbs(a1) < PetscAbs(a2) ? a1 : a2;

  *ox = a * px / pn;
  *oy = a * py / pn;
  
  if( a > sqrt(2) )
  {
    printf("a: %e  b: %e  c: %e  discr: %e  ox: %e  oy: %e\n", a, b, c, discr, *ox, *oy);
    printf("{");
    int x, y;
    for( y=0; y < 3; y++ )
      for( x=0; x < 3; x++ )
        printf("%f, ", phi[y][x]);
    printf("}\n");

    OrthoProjFirstOrder2D(phi,PHI_INF,ox,oy);
  }

  PetscFunctionReturn(0);
}
*/

#undef __FUNCT__
#define __FUNCT__ "OrthogonalProjection3D"
PetscInt EVENT_OrthogonalProjection3D;
PetscErrorCode OrthogonalProjection3D( double phi[3][3][3], Coor *op )
{
  int m;
  const int i = 1, j = 1, k = 1;
  double p, px, py, pz, px2, py2, pz2, pxx, pyy, pzz, pxy, pxz, pyz, pn;
  double c, a, f, df;
  
  PetscFunctionBegin;
  p  = phi[k][j][i];
  px = ( phi[k][j][i+1] - phi[k][j][i-1] ) / 2.;
  py = ( phi[k][j+1][i] - phi[k][j-1][i] ) / 2.;
  pz = ( phi[k+1][j][i] - phi[k-1][j][i] ) / 2.;
  px2 = px * px;
  py2 = py * py;
  pz2 = pz * pz;
  pxx = ( phi[k][j][i+1] - 2.*p + phi[k][j][i-1] );
  pyy = ( phi[k][j+1][i] - 2.*p + phi[k][j-1][i] );
  pzz = ( phi[k+1][j][i] - 2.*p + phi[k-1][j][i] );
  pxy = ( phi[k][j-1][i-1] + phi[k][j+1][i+1] - phi[k][j+1][i-1] - phi[k][j-1][i+1] ) / 4.;
  pxz = ( phi[k-1][j][i-1] + phi[k+1][j][i+1] - phi[k+1][j][i-1] - phi[k-1][j][i+1] ) / 4.; 
  pyz = ( phi[k-1][j-1][i] + phi[k+1][j+1][i] - phi[k+1][j-1][i] - phi[k-1][j+1][i] ) / 4.;
  pn = sqrt( px2 + py2 + pz2 );
  
  c = ( px2 * pxx +
        py2 * pyy +
        pz2 * pzz +
        2 * px * py * pxy +
        2 * px * pz * pxz +
        2 * py * pz * pyz ) / ( px2 + py2 + pz2 );

  // Newton's iteration
  a = 0;
  for ( m = 0; m < 5; ++m) {
    f = p + pn * a + 0.5*a*a*c;
    df = pn + a*c;
    a = a - f / df;
  }

  if( PetscAbs(a) < 1 ) {
    op->x = a * px / pn;
    op->y = a * py / pn;
    op->z = a * pz / pn;
  } else {
    OrthogonalProjection2D_1st( phi, PHI_INF, op );
  }

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "OrthoProjFirstOrder2D"
PetscErrorCode OrthoProjFirstOrder2D( double phi[3][3], double PHI_INF, double *ox, double *oy )
{
  const int nei[5][2] = {{1,0},{0,1},{-1,0},{0,-1},{1,0}};
  int i = 1, j = 1;
  int i1, j1, i2, j2;
  double p0, p1, p2, d1, d2, dd, b, min_dist;
  int k;
  
  PetscFunctionBegin;

  min_dist = PHI_INF;
  
  for( k = 0; k < 4; k++)
  {
    i1 = i + nei[k  ][0];
    j1 = j + nei[k  ][1];
    i2 = i + nei[k+1][0];
    j2 = j + nei[k+1][1];
    p0 = phi[j ][i ];
    p1 = phi[j1][i1];
    p2 = phi[j2][i2]; 

    if( p0 * p1 < 0. )
    {
      d1 = p0 / (p0 - p1);
      if( d1 < min_dist )
      {
        min_dist = d1;
        *ox = nei[k][0] * d1;
        *oy = nei[k][1] * d1;
      }
      if( p0 * p2 < 0. )
      {
        d2 = p0 / (p0 - p2);
        b = d1*d2 / (d1*d1 + d2*d2);
        dd = sqrt( (b*d1)*(b*d1) + (b*d2)*(b*d2) );
        if( dd < min_dist )
        {
          min_dist = dd;
          *ox = b * d1;
          *oy = b * d2;
        }
      }
    }
  }
  
  PetscFunctionReturn(0);
}
