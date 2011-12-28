#include "FiberField.h"

PetscErrorCode SphericalDistribution(PetscRandom rnd, Coor *d);
PetscErrorCode Rotation(Coor x0, Coor x1, Coor *n, Coor *r, Coor *s);
PetscErrorCode MatVec( Coor *n, Coor *r, Coor *s, Coor *a, Coor *d);

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc, char **args) {

  PetscErrorCode ierr;
  ierr = DCellInit(); CHKERRQ(ierr);

  FiberField fibers;
  ierr = FiberFieldCreate(PETSC_COMM_WORLD, &fibers); CHKERRQ(ierr);

  Vertex v1=0, v2;
  ierr = VertexCreate(fibers, &v1); CHKERRQ(ierr);

  PetscRandom rnd;
  ierr = PetscRandomCreate(PETSC_COMM_WORLD,&rnd); CHKERRQ(ierr);
  ierr = PetscRandomSetType(rnd,PETSCRAND48); CHKERRQ(ierr);

  Coor x0;
  Coor n,r,s,a,d;
  Edge e;
  int i;
  for (i = 0; i < 1e1; ++i) {
    ierr = VertexCreate(fibers, &v2); CHKERRQ(ierr);

    ierr = MemCacheAlloc(fibers->mcEdges,&e); CHKERRQ(ierr);
    ierr = UniqueIDGenerate(fibers->eid,&e->id); CHKERRQ(ierr);
    e->a = v1;
    e->id_a = v1->id;
    e->b = v2;
    e->id_b = v2->id;
    ierr = VertexAddEdge(v1,e); CHKERRQ(ierr);
    ierr = VertexAddEdge(v2,e); CHKERRQ(ierr);

    SphericalDistribution(rnd, &a);
    Rotation(x0,v1->X,&n,&r,&s);
    MatVec(&n,&r,&s,&a, &d);

    v2->X.x = v1->X.x + d.x;
    v2->X.y = v1->X.y + d.y;
    v2->X.z = v1->X.z + d.z;

    x0 = v1->X;
    v1 = v2;
  }

  ierr = FiberFieldDestroy(fibers); CHKERRQ(ierr);
	ierr = DCellFinalize(); CHKERRQ(ierr);
  return 0;
}

PetscErrorCode SphericalDistribution(PetscRandom rnd, Coor *d)
{
  PetscReal u,q;
  const PetscReal a = 0.95;

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
