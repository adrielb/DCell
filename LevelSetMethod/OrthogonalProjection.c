#include "LevelSetMethod.h"
#include "LSM_private.h"

#undef __FUNCT__
#define __FUNCT__ "OrthogonalProjection2D"
PetscErrorCode OrthogonalProjection2D( double phi3[3][3], double phi[5][5], Coor *op )
{
  Coor opB, opL;
  PetscReal distB, distL;
//  Coor opQ; 
//  PetscReal distQ;
  PetscFunctionBegin;
  OrthogonalProjection2D_Bicubic( phi, &opB, &distB );
//  OrthogonalProjection2D_Quadratic( phi3, &opQ, &distQ );
  OrthogonalProjection2D_Linear( phi3, &opL, &distL );

  if( distB < distL ) {
    *op = opB;
  } else {
    *op = opL;
  }
  PetscFunctionReturn(0);
}

// Some improvements of the fast marching method
// SIAM J Sci Comput Vol. 23, No. 1, pp 230-244

#undef __FUNCT__
#define __FUNCT__ "OrthogonalProjection2D_Bicubic"
PetscErrorCode OrthogonalProjection2D_Bicubic( double phi[5][5], Coor *op, PetscReal *minDist )
{
  const PetscReal TOL = 1e-3;
  const PetscReal TOL2 = TOL*TOL; // TOL^2
  const int MAXITER = 20;
  PetscReal a[16];
  PetscReal p, px, py, pp, xp;
  PetscReal x, x2, x3;
  PetscReal y, y2, y3;
  PetscReal norm;
  PetscReal dist, mindist = PETSC_MAX_REAL;
  Coor d1,d2;
  Coor x0;
  int i,j,m;

  for (j = 0; j <= 1; ++j) {
    for (i = 0; i <= 1; ++i) {
      // Bicubic interpolation
      a[0] = phi[1 + j][1 + i];
      a[1] = (-phi[1 + j][i] + phi[1 + j][2 + i])/2.;
      a[2] = phi[1 + j][i] - 3*phi[1 + j][1 + i] + 2*phi[1 + j][2 + i] + (phi[1 + j][1 + i] - phi[1 + j][3 + i])/2.;
      a[3] = 2*phi[1 + j][1 + i] - 2*phi[1 + j][2 + i] + (-phi[1 + j][i] + phi[1 + j][2 + i])/2. + (-phi[1 + j][1 + i] + phi[1 + j][3 + i])/2.;
      a[4] = (-phi[j][1 + i] + phi[2 + j][1 + i])/2.;
      a[5] = (phi[j][i] - phi[j][2 + i] - phi[2 + j][i] + phi[2 + j][2 + i])/4.;
      a[6] = (-3*(-phi[j][1 + i] + phi[2 + j][1 + i]))/2. + (-phi[j][i] + phi[j][2 + i] + phi[2 + j][i] - phi[2 + j][2 + i])/2. + (3*(-phi[j][2 + i] + phi[2 + j][2 + i]))/2. + (-phi[j][1 + i] + phi[j][3 + i] + phi[2 + j][1 + i] - phi[2 + j][3 + i])/4.;
      a[7] = -phi[j][1 + i] + phi[j][2 + i] + phi[2 + j][1 + i] - phi[2 + j][2 + i] + (phi[j][i] - phi[j][2 + i] - phi[2 + j][i] + phi[2 + j][2 + i])/4. + (phi[j][1 + i] - phi[j][3 + i] - phi[2 + j][1 + i] + phi[2 + j][3 + i])/4.;
      a[8] = phi[j][1 + i] - 3*phi[1 + j][1 + i] + 2*phi[2 + j][1 + i] + (phi[1 + j][1 + i] - phi[3 + j][1 + i])/2.;
      a[9] = (-3*(-phi[1 + j][i] + phi[1 + j][2 + i]))/2. + (-phi[j][i] + phi[j][2 + i] + phi[2 + j][i] - phi[2 + j][2 + i])/2. + (3*(-phi[2 + j][i] + phi[2 + j][2 + i]))/2. + (-phi[1 + j][i] + phi[1 + j][2 + i] + phi[3 + j][i] - phi[3 + j][2 + i])/4.;
      a[10] = phi[j][i] - phi[j][2 + i] + 9*phi[1 + j][1 + i] - 9*phi[1 + j][2 + i] + 3*(-phi[1 + j][i] + phi[1 + j][2 + i]) + (3*(-phi[1 + j][1 + i] + phi[1 + j][3 + i]))/2. - phi[2 + j][i] - 9*phi[2 + j][1 + i] + 3*(-phi[j][1 + i] + phi[2 + j][1 + i]) + 10*phi[2 + j][2 + i] - 3*(-phi[j][2 + i] + phi[2 + j][2 + i]) - 3*(-phi[2 + j][i] + phi[2 + j][2 + i]) - (3*(-phi[2 + j][1 + i] + phi[2 + j][3 + i]))/2. + (phi[j][1 + i] - phi[j][3 + i] - phi[2 + j][1 + i] + phi[2 + j][3 + i])/2. + (3*(-phi[1 + j][1 + i] + phi[3 + j][1 + i]))/2. - (3*(-phi[1 + j][2 + i] + phi[3 + j][2 + i]))/2. + (phi[1 + j][i] - phi[1 + j][2 + i] - phi[3 + j][i] + phi[3 + j][2 + i])/2. + (phi[1 + j][1 + i] - phi[1 + j][3 + i] - phi[3 + j][1 + i] + phi[3 + j][3 + i])/4.;
      a[11] = -5*phi[1 + j][1 + i] + 5*phi[1 + j][2 + i] - (3*(-phi[1 + j][i] + phi[1 + j][2 + i]))/2. - (3*(-phi[1 + j][1 + i] + phi[1 + j][3 + i]))/2. + 6*phi[2 + j][1 + i] - 2*(-phi[j][1 + i] + phi[2 + j][1 + i]) + (-phi[j][i] + phi[j][2 + i] + phi[2 + j][i] - phi[2 + j][2 + i])/2. - 6*phi[2 + j][2 + i] + 2*(-phi[j][2 + i] + phi[2 + j][2 + i]) + (3*(-phi[2 + j][i] + phi[2 + j][2 + i]))/2. + (-phi[j][1 + i] + phi[j][3 + i] + phi[2 + j][1 + i] - phi[2 + j][3 + i])/2. + (3*(-phi[2 + j][1 + i] + phi[2 + j][3 + i]))/2. - phi[3 + j][1 + i] + (-phi[1 + j][i] + phi[1 + j][2 + i] + phi[3 + j][i] - phi[3 + j][2 + i])/4. + phi[3 + j][2 + i] + (-phi[1 + j][1 + i] + phi[1 + j][3 + i] + phi[3 + j][1 + i] - phi[3 + j][3 + i])/4.;
      a[12] = 2*phi[1 + j][1 + i] - 2*phi[2 + j][1 + i] + (-phi[j][1 + i] + phi[2 + j][1 + i])/2. + (-phi[1 + j][1 + i] + phi[3 + j][1 + i])/2.;
      a[13] = -phi[1 + j][i] + phi[1 + j][2 + i] + phi[2 + j][i] - phi[2 + j][2 + i] + (phi[j][i] - phi[j][2 + i] - phi[2 + j][i] + phi[2 + j][2 + i])/4. + (phi[1 + j][i] - phi[1 + j][2 + i] - phi[3 + j][i] + phi[3 + j][2 + i])/4.;
      a[14] = -5*phi[1 + j][1 + i] + 6*phi[1 + j][2 + i] - 2*(-phi[1 + j][i] + phi[1 + j][2 + i]) - phi[1 + j][3 + i] + 5*phi[2 + j][1 + i] - (3*(-phi[j][1 + i] + phi[2 + j][1 + i]))/2. + (-phi[j][i] + phi[j][2 + i] + phi[2 + j][i] - phi[2 + j][2 + i])/2. - 6*phi[2 + j][2 + i] + (3*(-phi[j][2 + i] + phi[2 + j][2 + i]))/2. + 2*(-phi[2 + j][i] + phi[2 + j][2 + i]) + (-phi[j][1 + i] + phi[j][3 + i] + phi[2 + j][1 + i] - phi[2 + j][3 + i])/4. + phi[2 + j][3 + i] - (3*(-phi[1 + j][1 + i] + phi[3 + j][1 + i]))/2. + (-phi[1 + j][i] + phi[1 + j][2 + i] + phi[3 + j][i] - phi[3 + j][2 + i])/2. + (3*(-phi[1 + j][2 + i] + phi[3 + j][2 + i]))/2. + (-phi[1 + j][1 + i] + phi[1 + j][3 + i] + phi[3 + j][1 + i] - phi[3 + j][3 + i])/4.;
      a[15] = -phi[j][1 + i] + phi[j][2 + i] - phi[1 + j][i] + 2*phi[1 + j][1 + i] - 2*phi[1 + j][2 + i] + phi[1 + j][3 + i] + phi[2 + j][i] - 2*phi[2 + j][1 + i] + 2*phi[2 + j][2 + i] + (phi[j][i] - phi[j][2 + i] - phi[2 + j][i] + phi[2 + j][2 + i])/4. - phi[2 + j][3 + i] + (phi[j][1 + i] - phi[j][3 + i] - phi[2 + j][1 + i] + phi[2 + j][3 + i])/4. + phi[3 + j][1 + i] - phi[3 + j][2 + i] + (phi[1 + j][i] - phi[1 + j][2 + i] - phi[3 + j][i] + phi[3 + j][2 + i])/4. + (phi[1 + j][1 + i] - phi[1 + j][3 + i] - phi[3 + j][1 + i] + phi[3 + j][3 + i])/4.;

      // Newton's iteration
      x0.x = 1-i;
      x0.y = 1-j;
      x = x0.x;
      y = x0.y;
//printf("{{%f,%f}, ", x+i-1, y+j-1);
//printf("%e \n",norm);
      for (m = 0; m < MAXITER; ++m) {
        x2 = x*x;
        x3 = x2*x;
        y2 = y*y;
        y3 = y2*y;
        p = a[0] + x*a[1] + x2*a[2] + x3*a[3] + y*a[4] + x*y*a[5] + x2*y*a[6] + x3*y*a[7] + y2*a[8] + x*y2*a[9] + x2*y2*a[10] + x3*y2*a[11] + y3*a[12] + x*y3*a[13] + x2*y3*a[14] + x3*y3*a[15];
        px = a[1] + 2*x*a[2] + 3*x2*a[3] + y*a[5] + 2*x*y*a[6] + 3*x2*y*a[7] + y2*a[9] + 2*x*y2*a[10] + 3*x2*y2*a[11] + y3*a[13] + 2*x*y3*a[14] + 3*x2*y3*a[15];
        py = a[4] + x*a[5] + x2*a[6] + x3*a[7] + 2*y*a[8] + 2*x*y*a[9] + 2*x2*y*a[10] + 2*x3*y*a[11] + 3*y2*a[12] + 3*x*y2*a[13] + 3*x2*y2*a[14] + 3*x3*y2*a[15];

        // Solves p( y ) == 0
        pp = px*px + py*py;
        d1.x = -p * px / pp;
        d1.y = -p * py / pp;

        // Solves dp(y) . (x0 - y) == 0
        xp = (x0.x - x) * px +
             (x0.y - y) * py;
        d2.x = (x0.x - x) - xp / pp * px;
        d2.y = (x0.y - y) - xp / pp * py;

        x += d1.x + 0.5*d2.x;
        y += d1.y + 0.5*d2.y;
        norm = d1.x*d1.x + d1.y*d1.y +
               d2.x*d2.x + d2.y*d2.y ;
//printf("{%f,%f}, ", x+i-1, y+j-1);
//printf("%e \n",norm);
        if( norm < TOL2 ) break;
      }
//printf("},\n");
//printf("m = %d\n", m);
      if( m == MAXITER ) { // convergence failed
        dist = PETSC_MAX_REAL;
      } else {
        if( x < -TOL || 1+TOL < x ||
            y < -TOL || 1+TOL < y  ) { // converged outside of domain
          dist = PETSC_MAX_REAL;
        } else {
          dist = PetscSqrtScalar(
              (x - x0.x) * (x - x0.x) +
              (y - x0.y) * (y - x0.y) );
        }
      }
      if( dist < mindist ) {
        mindist = dist;
        op->x = x + i - 1;
        op->y = y + j - 1;
        op->z = 0;
      }
    } // i
  } // j

  *minDist = mindist;
  return 0;
}

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
#define __FUNCT__ "OrthogonalProjection2D_Quadratic"
PetscErrorCode OrthogonalProjection2D_Quadratic( double phi[3][3], Coor *op, PetscReal *mindist )
{
  int m;
  const int i = 1, j = 1; //makes indexing easier
  double p, px, py, px2, py2, pxx, pyy, pxy, pn, c, a, f, df;

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

  *mindist = 2;
  if( PetscAbs(a) < 1 ) {
    op->x = a * px / pn;
    op->y = a * py / pn;
    *mindist = a;
  }
  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "OrthogonalProjection2D_Linear"
PetscErrorCode OrthogonalProjection2D_Linear( double phi[3][3],  Coor *op, PetscReal *mindist )
{
  int i;
  int x,y;
  double dist;

  *mindist = 2;
  for (i = 0; i < 4; ++i) {
    x = STAR[i][0] + 1;
    y = STAR[i][1] + 1;
    if( phi[1][1] * phi[y][x] <= 0 ) {
      dist = phi[1][1] / (phi[1][1] - phi[y][x]);
      if( dist < *mindist ) {
        *mindist = dist;
        op->x = dist * STAR[i][0];
        op->y = dist * STAR[i][1];
      }
    }
  }
  return 0;
}
