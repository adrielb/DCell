#include "LevelSetMethod.h"
#include "LSM_private.h"

#undef __FUNCT__
#define __FUNCT__ "OrthogonalProjection3D"
PetscErrorCode OrthogonalProjection3D( double phi3[3][3][3], double phi[5][5][5], Coor *op ) {
  const PetscReal TOL = 1e-3;
  const PetscReal TOL2 = TOL*TOL; // TOL^2
  const int MAXITER = 20;
  Coor d1,d2;
  Coor x0;
  PetscReal x, x2, x3;
  PetscReal y, y2, y3;
  PetscReal z, z2, z3;
  PetscReal p, px, py, pz;
  PetscReal pp, xp;
  PetscReal norm;
  PetscReal dist, mindist = PETSC_MAX_REAL;
  int i, j, k, m;
  PetscReal a[64]; // 4x4x4 coeff in tricubic interpolation

  PetscFunctionBegin;
  for (k = 0; k <= 1; ++k) {
    for (j = 0; j <= 1; ++j) {
      for (i = 0; i <= 1; ++i) {

        #include "TricubicInverse.h"
        // a[0] = ...
        // ...
        // a[63]= ...

        // Newton's iteration
        x0.x = 1-i;
        x0.y = 1-j;
        x0.z = 1-k;
        x = x0.x;
        y = x0.y;
        z = x0.z;

        for (m = 0; m < MAXITER; ++m) {
          x2 = x*x;
          x3 = x2*x;
          y2 = y*y;
          y3 = y2*y;
          z2 = z*z;
          z3 = z2*z;

          #include "TricubicInterpolation.h"
          // p  = ...
          // px = ...
          // py = ...
          // pz = ...

          // Solves p( y ) == 0
          pp = px*px + py*py + pz*pz;
          d1.x = -p * px / pp;
          d1.y = -p * py / pp;
          d1.z = -p * pz / pp;

          // Solves dp(y) . (x0 - y) == 0
          xp = (x0.x - x) * px +
               (x0.y - y) * py +
               (x0.z - z) * pz;
          d2.x = (x0.x - x) - xp / pp * px;
          d2.y = (x0.y - y) - xp / pp * py;
          d2.z = (x0.z - z) - xp / pp * pz;

          x += d1.x + 0.5*d2.x;
          y += d1.y + 0.5*d2.y;
          z += d1.z + 0.5*d2.z;
          norm = d1.x*d1.x + d1.y*d1.y + d1.z*d1.z +
                 d2.x*d2.x + d2.y*d2.y + d2.z*d2.z;
  //printf("{%f,%f}, ", x+i-1, y+j-1);
  //printf("%e \n",norm);
          if( norm < TOL2 ) break;
        }
  //printf("},\n");
  //printf("m = %d\n", m);
        if( m == MAXITER ) { // convergence failed
          dist = PETSC_MAX_REAL;
        } else {
          // converged outside of domain
          if( x < -TOL || 1+TOL < x ||
              y < -TOL || 1+TOL < y ||
              z < -TOL || 1+TOL < z ) {
            dist = PETSC_MAX_REAL;
          } else {
            dist = PetscSqrtScalar(
                (x - x0.x) * (x - x0.x) +
                (y - x0.y) * (y - x0.y) +
                (z - x0.z) * (z - x0.z) );
          }
        }
        if( dist < mindist ) {
          mindist = dist;
          op->x = x + i - 1;
          op->y = y + j - 1;
          op->z = z + k - 1;
        }
      } // i
    } // j
  } // k
//  printf("MINDIST: %e\n", mindist);
//  printf("Local: {%d, %d}\n", i, j);
  if( PetscAbs( mindist - PETSC_MAX_REAL) < 1 ) {
    OrthogonalProjection3D_Quadratic( phi3, op);
/*
    printf("PHI\n");
    printf("{");
    for (j = 0; j < 5; ++j) {
      printf("{");
      for (i = 0; i < 4; ++i) {
        printf("%f, ", phi[j][i]);
      }
      printf("%f},", phi[j][4]);
    }
    printf("}\n");
    printf("{");
    for (i = 0; i < 15; ++i) {
      printf("%f, ", a[i]);
    }
    printf("%f}\n", a[15]);
    exit(1);
*/
  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "OrthogonalProjection3D_Quadratic"
PetscErrorCode OrthogonalProjection3D_Quadratic( double phi[3][3][3], Coor *op)
{
  int m;
  const int i = 1, j = 1, k = 1; //makes indexing easier
  double px,  py,  pz,
         px2, py2, pz2,
         pxx, pyy, pzz,
         pxy, pxz, pyz;
  double p, pn, c, a, f, df;

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
        2 * py * pz * pyz
      ) / ( px2+py2+pz2 );

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
    OrthogonalProjection3D_Linear( phi, op );
  }

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "OrthogonalProjection3D_Linear"
PetscErrorCode OrthogonalProjection3D_Linear( double phi[3][3][3],  Coor *op )
{
  int i, x, y, z;
  double mindist = 2, dist;

  PetscFunctionBegin;
  for (i = 0; i < 6; ++i) {
    x = STAR[i][0] + 1;
    y = STAR[i][1] + 1;
    z = STAR[i][2] + 1;
    if( phi[1][1][1] * phi[z][y][x] <= 0 ) {
      dist = phi[1][1][1] / (phi[1][1][1] - phi[z][y][x]);
      if( dist < mindist ) {
        mindist = dist;
        op->x = dist * STAR[i][0];
        op->y = dist * STAR[i][1];
        op->z = dist * STAR[i][2];
      }
    }
  }
  PetscFunctionReturn(0);
}
