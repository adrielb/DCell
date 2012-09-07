#include "ImmersedInterfaceMethod.h"
#include "IIM_private.h"

#undef __FUNCT__
#define __FUNCT__ "LocalCoorCreate"
PetscErrorCode LocalCoorCreate( LocalCoor *lc )
{
  PetscErrorCode ierr;
  LocalCoor l;
  
  PetscFunctionBegin;

  ierr = PetscNew( struct _LocalCoor, &l); CHKERRQ(ierr);
  ierr = ArrayCreate("localcoor", 3*sizeof(PetscReal), &l->coor); CHKERRQ(ierr);

  *lc = l;
  
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "LocalCoorDestroy"
PetscErrorCode LocalCoorDestroy( LocalCoor lc )
{
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  ierr = ArrayDestroy(lc->coor); CHKERRQ(ierr);
  ierr = PetscFree( lc ); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

void LocalCoorSetLength( LocalCoor lc, PetscInt len)
{
  ArraySetSize(lc->coor, len);
}

void LocalCoor3DTangential( IIMIrregularNode *n )
{
  PetscReal h;

  if( n->nx * n->nx + n->ny * n->ny >=
      n->nx * n->nx + n->nz * n->nz )
  {
    h = sqrt( n->nx * n->nx + n->ny * n->ny );
    n->sx =  n->ny / h;
    n->sy = -n->nx / h;
    n->sz = 0;

    n->rx =  n->nx * n->nz;
    n->ry =  n->ny * n->nz;
    n->rz = -n->nx * n->nx - n->ny * n->ny;
  } else {
    h = sqrt( n->nx * n->nx + n->nz * n->nz );
    n->sx = n->nz / h;
    n->sy = 0;
    n->sz =-n->nx / h;

    n->rx = -n->nx * n->ny;
    n->ry =  n->nx * n->nx + n->nz * n->nz;
    n->rz = -n->ny * n->nz;
  }
  h = sqrt( n->rx * n->rx + n->ry * n->ry + n->rz * n->rz );
  n->rx /= h;
  n->ry /= h;
  n->rz /= h;
}

void LocalCoorSolve( LocalCoor lc, Coor dh, IIMIrregularNode *N, IIMIrregularNode *nodes[] )
{
  int i;
  PetscReal x, y, z;
  PetscReal *s, *r, *n;
  LocalCoorGetVecs( lc, &s, &n, &r);
  const int len = ArrayLength(lc->coor);

  for ( i = 0; i < len; ++i)
  {
    x = (nodes[i]->X.x - N->X.x) * dh.x;   // s <--> x
    y = (nodes[i]->X.y - N->X.y) * dh.y;   // n <--> y
    z = (nodes[i]->X.z - N->X.z) * dh.z;   // r <--> z
    s[i] = N->sx * x + N->sy * y + N->sz * z;
    n[i] = N->nx * x + N->ny * y + N->nz * z;
    r[i] = N->rx * x + N->ry * y + N->rz * z;

    s[i] = PetscSign(s[i]) * PetscSqrtReal( s[i]*s[i] + n[i]*n[i] );
    r[i] = PetscSign(r[i]) * PetscSqrtReal( r[i]*r[i] + n[i]*n[i] );
  }
}

//                                               x              y              z
void LocalCoorGetVecs( LocalCoor lc, PetscReal **s, PetscReal **n, PetscReal **r )
{
  const int len = ArrayLength(lc->coor);
  PetscReal *coor = ArrayGetData(lc->coor);
  if( s )  *s = &coor[0*len];
  if( n )  *n = &coor[1*len];
  if( r )  *r = &coor[2*len];
}
