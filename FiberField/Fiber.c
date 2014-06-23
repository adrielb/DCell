#include "FiberField.h"
#include "FiberField_private.h"

PetscErrorCode BoundaryCheck( Coor min, Coor max, Coor *X, PetscBool *hitBoundary );
PetscErrorCode SphericalDistribution(PetscRandom rnd, PetscReal a, Coor *d);
PetscErrorCode Rotation(Coor x0, Coor x1, Coor *n, Coor *r, Coor *s);
PetscErrorCode MatVec( Coor *n, Coor *r, Coor *s, Coor *a, Coor *d);

#undef __FUNCT__
#define __FUNCT__ "FiberFieldInitLocalFiber"
PetscErrorCode FiberFieldInitLocalFiber( FiberField fibers, int numVerticies, PetscReal l0,
    FiberTypeID vertexType, FiberTypeID edgeType, FiberTypeID bendingEdgeType )
{
  int i;
  int flen = numVerticies;
  Array fiber = fibers->fiber;
  Vertex v0;
  Vertex v1;
  Vertex v2;
  PetscRandom rnd = fibers->rnd;
  PetscBool hitBoundary;
  Coor rndSphere;
  Coor n;
  Coor r;
  Coor s;
  Coor d;
  Coor lmin = fibers->localBounds.min;
  Coor lmax = fibers->localBounds.max;
  PetscReal NORTH_HEMISPHERE = 0.65;
  const PetscReal WHOLE_SPHERE = 0.0;
  const PetscReal DELTA = PETSC_SMALL;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = ArraySetSize(fiber, 0); CHKERRQ(ierr);
  ierr = ArraySetMaxSize( fibers->verts, flen + ArrayLength(fibers->verts) ); CHKERRQ(ierr);

  // add small delta
  lmin.x += DELTA;
  lmin.y += DELTA;
  lmin.z += DELTA;
  lmax.x -= DELTA;
  lmax.y -= DELTA;
  lmax.z -= DELTA;

  ierr = PetscOptionsGetReal(0,"-fiber_bend",&NORTH_HEMISPHERE,0); CHKERRQ(ierr);

  // Initialize v0 anywhere in the local bbox
  ierr = FiberFieldAddVertex( fibers, vertexType, &v0); CHKERRQ(ierr);
  ierr = ArrayAppendPtr( fiber, v0); CHKERRQ(ierr);
  PetscRandomGetValue(rnd,&v0->X.x); // [0,1]
  PetscRandomGetValue(rnd,&v0->X.y); // [0,1]
  PetscRandomGetValue(rnd,&v0->X.z); // [0,1]
  // map [0, 1] to local bbox [lmin, lmax]
  v0->X.x = (1-v0->X.x) * lmin.x + v0->X.x * lmax.x;
  v0->X.y = (1-v0->X.y) * lmin.y + v0->X.y * lmax.y;
  v0->X.z = (1-v0->X.z) * lmin.z + v0->X.z * lmax.z;

  // Initialize v1 anywhere spherically from v0
  ierr = FiberFieldAddVertex( fibers, vertexType, &v1); CHKERRQ(ierr);
  ierr = ArrayAppendPtr( fiber, v1); CHKERRQ(ierr);
  SphericalDistribution(rnd, WHOLE_SPHERE, &d);
  v1->X.x = v0->X.x + l0*d.x;
  v1->X.y = v0->X.y + l0*d.y;
  v1->X.z = v0->X.z + l0*d.z;
  BoundaryCheck( lmin, lmax, &v1->X, &hitBoundary );
  if (hitBoundary) {
    ierr = FiberFieldAddEdge( fibers, v0, v1, edgeType, l0 ); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

  for (i = 1; i < flen-1; i++ ) {
    ierr = FiberFieldAddVertex( fibers, vertexType, &v2); CHKERRQ(ierr);
    ierr = ArrayAppendPtr( fiber, v2); CHKERRQ(ierr);

    SphericalDistribution(rnd, NORTH_HEMISPHERE, &rndSphere);
    Rotation(v0->X,v1->X,&n,&r,&s);
    MatVec(&n,&r,&s,&rndSphere, &d);

    v2->X.x = v1->X.x + l0*d.x;
    v2->X.y = v1->X.y + l0*d.y;
    v2->X.z = v1->X.z + l0*d.z;

    BoundaryCheck( lmin, lmax, &v2->X, &hitBoundary );
    if (hitBoundary) {
      break;
    }

    v0 = v1;
    v1 = v2;
  }

  flen = ArrayLength( fiber );
  Vertex *v = ArrayGetData( fiber );
  // Link fiber edges: v_i -- v_i+1
  for (i = 0; i < flen-1; i++) {
    ierr = FiberFieldAddEdge( fibers, v[i], v[i+1], edgeType, l0 ); CHKERRQ(ierr);
  }

  // Link bending edges: v_i-1 -- v_i+1
  for (i = 1; i < flen-1; i++) {
    ierr = FiberFieldAddEdge( fibers, v[i-1], v[i+1], bendingEdgeType, 2*l0 ); CHKERRQ(ierr);
  }

  PetscFunctionReturn(0);
}

// http://mathworld.wolfram.com/SpherePointPicking.html
PetscErrorCode SphericalDistribution(PetscRandom rnd, PetscReal a, Coor *d)
{
  PetscReal u,q;

  PetscRandomGetValue(rnd,&u); // [0,1]
  u = (1-a) * u + a; // [a,1]
  PetscRandomGetValue(rnd,&q); // [0,1]
  q = 2 * PETSC_PI * q; // [0, 2pi]

  d->x = PetscSqrtScalar( 1 - u*u ) * cos(q);
  d->y = PetscSqrtScalar( 1 - u*u ) * sin(q);
  d->z = u;

  return 0;
}

PetscErrorCode Rotation(Coor x0, Coor x1, Coor *n, Coor *r, Coor *s)
{
  PetscReal mag;
  n->x = x1.x - x0.x;
  n->y = x1.y - x0.y;
  n->z = x1.z - x0.z;
  mag = PetscSqrtScalar(n->x*n->x + n->y*n->y + n->z*n->z);
  n->x /= mag;
  n->y /= mag;
  n->z /= mag;

  PetscReal nx2 = n->x*n->x,
            ny2 = n->y*n->y,
            nz2 = n->z*n->z;
  if( nx2+ ny2 > nx2 + nz2 ) {
    mag = PetscSqrtScalar(nx2 + ny2);
    r->x =  n->y / mag;
    r->y = -n->x / mag;
    r->z = 0;

    s->x = n->x * n->z;
    s->y = n->y * n->z;
    s->z = -nx2-ny2;
  } else {
    mag = PetscSqrtScalar(nx2 + nz2);
    r->x =  n->z / mag;
    r->y = 0;
    r->z = -n->x / mag;

    s->x = -n->x * n->y;
    s->y = nx2 + nz2;
    s->z = -n->y * n->z;
  }
  mag = PetscSqrtScalar(s->x*s->x + s->y*s->y + s->z*s->z);
  s->x /= mag;
  s->y /= mag;
  s->z /= mag;

  return 0;
}

PetscErrorCode MatVec( Coor *n, Coor *r, Coor *s, Coor *a, Coor *d)
{
  d->x = a->z*n->x + a->x*r->x + a->y*s->x;
  d->y = a->z*n->y + a->x*r->y + a->y*s->y;
  d->z = a->z*n->z + a->x*r->z + a->y*s->z;
  return 0;
}

PetscErrorCode BoundaryCheck( Coor min, Coor max, Coor *X, PetscBool *hitBoundary )
{ 
  *hitBoundary = PETSC_FALSE;

  // upper boundary 
  if (X->x > max.x) { X->x = max.x; *hitBoundary = PETSC_TRUE; }
  if (X->y > max.y) { X->y = max.y; *hitBoundary = PETSC_TRUE; }
  if (X->z > max.z) { X->z = max.z; *hitBoundary = PETSC_TRUE; }

  // lower boundary
  if (X->x < min.x) { X->x = min.x; *hitBoundary = PETSC_TRUE; }
  if (X->y < min.y) { X->y = min.y; *hitBoundary = PETSC_TRUE; }
  if (X->z < min.z) { X->z = min.z; *hitBoundary = PETSC_TRUE; }

  return 0;
}
