#include "FiberField"

// p is your point, p.x is the x coord, p.y is the y coord
if (p.x < min.x || p.x > max.x ||
    p.y < min.y || p.y > max.y ||
    p.z < min.z || p.z > max.z ) {
    // Definitely not within the polygon!
}

inside = PETSC_FALSE;
for (k = 0; k < 2; ++k) {
  for (j = 0; j < 2; ++j) {
    for (i = 0; i < 2; ++i) {
      p.x = i == 0 ? min.x : max.x;
      p.y = k == 0 ? min.y : max.y;
      p.z = j == 0 ? min.z : max.z;
      inside = inBBox( p, bbox );
      if( inside )
        break;
    }
  }
}


FiberField_CollisionDetection( FiberField f )
{
  Coor X1, X2, X3, X4;
  Coor V1, V2, V3, V4;

  // xij = xi - xj
  Coor x21 = { X2.x - X1.x, X2.y - X1.y, X2.z - X1.z};
  Coor x31 = { X3.x - X1.x, X3.y - X1.y, X3.z - X1.z};
  Coor x41 = { X4.x - X1.x, X4.y - X1.y, X4.z - X1.z};
  Coor x43 = { X4.x - X3.x, X4.y - X3.y, X4.z - X3.z};
  Coor v21 = { X2.x - X1.x, X2.y - X1.y, X2.z - X1.z};
  Coor v31 = { X3.x - X1.x, X3.y - X1.y, X3.z - X1.z};
  Coor v41 = { X4.x - X1.x, X4.y - X1.y, X4.z - X1.z};

  // parallel edges
  // | x21 - x43 | < tol
  a = -(x21.z*x43.y) + x21.y*x43.z;
  b =   x21.z*x43.x - x21.x*x43.z;
  c = -(x21.y*x43.x) + x21.x*x43.y;
  if( sqrt(a*a+b*b+c*c) < tol )
    // is parallel
  else

  // a + b*t + c*t^2 + d*t^3 == 0
  a = -(x21.z*x31.y*x41.x) + x21.y*x31.z*x41.x + x21.z*x31.x*x41.y - x21.x*x31.z*x41.y
      - x21.y*x31.x*x41.z + x21.x*x31.y*x41.z;
  b = -(v41.z*x21.y*x31.x) + v41.y*x21.z*x31.x + v41.z*x21.x*x31.y - v41.x*x21.z*x31.y
      - v41.y*x21.x*x31.z + v41.x*x21.y*x31.z + v31.z*x21.y*x41.x
      - v31.y*x21.z*x41.x - v21.z*x31.y*x41.x + v21.y*x31.z*x41.x - v31.z*x21.x*x41.y
      + v31.x*x21.z*x41.y + v21.z*x31.x*x41.y - v21.x*x31.z*x41.y
      + v31.y*x21.x*x41.z - v31.x*x21.y*x41.z - v21.y*x31.x*x41.z + v21.x*x31.y*x41.z;
  c = -(v31.z*v41.y*x21.x) + v31.y*v41.z*x21.x + v31.z*v41.x*x21.y - v31.x*v41.z*x21.y
      - v31.y*v41.x*x21.z + v31.x*v41.y*x21.z + v21.z*v41.y*x31.x
      - v21.y*v41.z*x31.x - v21.z*v41.x*x31.y + v21.x*v41.z*x31.y + v21.y*v41.x*x31.z
      - v21.x*v41.y*x31.z - v21.z*v31.y*x41.x + v21.y*v31.z*x41.x
      + v21.z*v31.x*x41.y - v21.x*v31.z*x41.y - v21.y*v31.x*x41.z + v21.x*v31.y*x41.z;
  d = -(v21.z*v31.y*v41.x) + v21.y*v31.z*v41.x + v21.z*v31.x*v41.y - v21.x*v31.z*v41.y
      - v21.y*v31.x*v41.z + v21.x*v31.y*v41.z;

  // first root in [0, dt]
  int numIC = 10;
  int numNewton = 10;
  for (i = 0; i < numIC; ++i) {
    t = (dt * i) / ( numIC - 1. );
    for (m = 0; m < numNewton; ++m) {
      f = a + b*t + c*t*t + d*t*t*t;
      df = b + 2*c*t + 3*d*t*t;
    }
  }

  if( f )

  Ic = m*vn / 2;
  Ir = -PetscMin(dt*k*d, mass*(0.1*d/dt - vn));
  It = 2*I / ( a*a + (1-a)*(1-a) + b*b + (1-b)*(1-b) );
  V1.x = V1.x + (1-a)(It/m)
}
