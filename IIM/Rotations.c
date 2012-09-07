#include "ImmersedInterfaceMethod.h"
#include "IIM_private.h"
//TODO: eliminate cross derivatives from rotation application, dont need them for correction terms

/*
GlobalToLocal
{
  x * nx + y * ny + z * nz;
  x * sx + y * sy + z * sz;
  x * rx + y * ry + z * rz;
}
*/

void IIMRotation1st_2D( IIMIrregularNode *n, Jump *j )
{
  j->x = j->e*n->nx + j->n*n->sx;
  j->y = j->e*n->ny + j->n*n->sy;
}

// Local to global coordinate rotation for 1st order terms
void IIMLocalToGlobal_1st( IIMIrregularNode *n, Jump *j)
{
  j->x = j->e*n->nx + j->t*n->rx + j->n*n->sx;
  j->y = j->e*n->ny + j->t*n->ry + j->n*n->sy;
  j->z = j->e*n->nz + j->t*n->rz + j->n*n->sz;
}

void IIMRotation2nd_2D( IIMIrregularNode *n, Jump *j )
{
  j->xx = n->nx*(j->ee*n->nx + j->ne*n->sx) + n->sx*(j->ne*n->nx + j->nn*n->sx);
  j->xy = n->ny*(j->ee*n->nx + j->ne*n->sx) + n->sy*(j->ne*n->nx + j->nn*n->sx);
  j->yy = n->ny*(j->ee*n->ny + j->ne*n->sy) + n->sy*(j->ne*n->ny + j->nn*n->sy);
}

// Local to global coordinate rotation for 2nd order terms
void IIMLocalToGlobal_2nd( IIMIrregularNode *n, Jump *j)
{
  j->xx = n->nx*(j->ee*n->nx + j->te*n->rx + j->ne*n->sx) + n->sx*(j->ne*n->nx + j->nt*n->rx + j->nn*n->sx) + n->rx*(j->te*n->nx + j->tt*n->rx + j->nt*n->sx);
  j->xy = n->ny*(j->ee*n->nx + j->te*n->rx + j->ne*n->sx) + n->ry*(j->te*n->nx + j->tt*n->rx + j->nt*n->sx) + (j->ne*n->nx + j->nt*n->rx + j->nn*n->sx)*n->sy;
  j->xz = n->nz*(j->ee*n->nx + j->te*n->rx + j->ne*n->sx) + n->rz*(j->te*n->nx + j->tt*n->rx + j->nt*n->sx) + (j->ne*n->nx + j->nt*n->rx + j->nn*n->sx)*n->sz;
  j->yy = n->ny*(j->ee*n->ny + j->te*n->ry + j->ne*n->sy) + n->sy*(j->ne*n->ny + j->nt*n->ry + j->nn*n->sy) + n->ry*(j->te*n->ny + j->tt*n->ry + j->nt*n->sy);
  j->yz = n->nz*(j->ee*n->ny + j->te*n->ry + j->ne*n->sy) + n->rz*(j->te*n->ny + j->tt*n->ry + j->nt*n->sy) + (j->ne*n->ny + j->nt*n->ry + j->nn*n->sy)*n->sz;
  j->zz = n->nz*(j->ee*n->nz + j->te*n->rz + j->ne*n->sz) + n->sz*(j->ne*n->nz + j->nt*n->rz + j->nn*n->sz) + n->rz*(j->te*n->nz + j->tt*n->rz + j->nt*n->sz);
}
void IIMGlobalToLocal_1st( IIMIrregularNode *n, Jump *j)
{
  j->e = j->x*n->nx + j->y*n->ny + j->z*n->nz;
  j->n = j->x*n->sx + j->y*n->sy + j->z*n->sz;
  j->t = j->x*n->rx + j->y*n->ry + j->z*n->rz;
}
void IIMGlobalToLocal_2nd( IIMIrregularNode *n, Jump *j)
{
  j->ee = n->nx*(j->xx*n->nx + j->xy*n->ny + j->xz*n->nz) + n->ny*(j->xy*n->nx + j->yy*n->ny + j->yz*n->nz) + n->nz*(j->xz*n->nx + j->yz*n->ny + j->zz*n->nz);
  j->ne = (j->xx*n->nx + j->xy*n->ny + j->xz*n->nz)*n->sx + (j->xy*n->nx + j->yy*n->ny + j->yz*n->nz)*n->sy + (j->xz*n->nx + j->yz*n->ny + j->zz*n->nz)*n->sz;
  j->te = (j->xx*n->nx + j->xy*n->ny + j->xz*n->nz)*n->rx + (j->xy*n->nx + j->yy*n->ny + j->yz*n->nz)*n->ry + (j->xz*n->nx + j->yz*n->ny + j->zz*n->nz)*n->rz;
  j->nn = n->sx*(j->xx*n->sx + j->xy*n->sy + j->xz*n->sz) + n->sy*(j->xy*n->sx + j->yy*n->sy + j->yz*n->sz) + n->sz*(j->xz*n->sx + j->yz*n->sy + j->zz*n->sz);
  j->nt = n->rx*(j->xx*n->sx + j->xy*n->sy + j->xz*n->sz) + n->ry*(j->xy*n->sx + j->yy*n->sy + j->yz*n->sz) + n->rz*(j->xz*n->sx + j->yz*n->sy + j->zz*n->sz);
  j->tt = n->rx*(j->xx*n->rx + j->xy*n->ry + j->xz*n->rz) + n->ry*(j->xy*n->rx + j->yy*n->ry + j->yz*n->rz) + n->rz*(j->xz*n->rx + j->yz*n->ry + j->zz*n->rz);
}
