#include "ImmersedInterfaceMethod.h"

typedef struct
{
  PetscReal d, k, jp, jpn, jpss, jf;
} Jump;

void JumpConditionPressure2D( IrregularNode *n, Jump *j)
{
  j->d   = sqrt( PetscSqr( n->ox ) + PetscSqr( n->oy ) );
  j->k   = n->k;
  j->jp  = n->f1;
  j->jpn = n->df2;
  j->jpss= n->ddf1;
  j->jf  = 0;
}

void JumpConditionVelocity2D_X( IrregularNode *n, Jump *j)
{
  PetscReal mu = 1; //TODO: Temporary, until fluid property encapsulation worked out
  j->d    = sqrt( PetscSqr( n->ox ) + PetscSqr( n->oy ) );
  j->k    = n->k;
  j->jp   = 0;
  j->jpn  = n->f2 * n->ny / mu;
  j->jpss = 0;
  j->jf   = ( n->ny * n->df2 + n->nx * n->df1 ) / mu;
}

void JumpConditionVelocity2D_Y( IrregularNode *n, Jump *j)
{
  PetscReal mu = 1; //TODO: Temporary, until fluid property encapsulation worked out
  j->d    = sqrt( PetscSqr( n->ox ) + PetscSqr( n->oy ) );
  j->k    = n->k;
  j->jp   = 0;
  j->jpn  = -n->f2 * n->nx / mu;
  j->jpss = 0;
  j->jf   = ( n->nx * n->df2 - n->ny * n->df1 ) / mu;
}

void JumpConditionPressure3D( IrregularNode *n, Jump *j)
{
  j->d   = sqrt( PetscSqr( n->ox ) + PetscSqr( n->oy ) + PetscSqr( n->oz ) );
  j->k   = n->k;
  j->jp  = n->f1;
//  j->jpn = n->dFtds + n->dFtdr;
//  j->jpss= n->ddFndss + n->ddFndrr; // TODO: include second tangential direction in IrregularNodes
  j->jf  = 0;
  SETERRQ(1,"Not Implemented");
}

void JumpConditionVelocity3D_X( IrregularNode *n, Jump *j)
{
  j->d    = sqrt( PetscSqr( n->ox ) + PetscSqr( n->oy ) + PetscSqr( n->oz ) );
  j->k    = n->k;
  j->jp   = 0;
  j->jpn  = -n->f2 * n->sx;
  j->jpss = 0;
  /*
  j->jf   = j->jps*n->nz*n->ry - j->jps*n->ny*n->rz - j->jpr*n->nz*n->sy + 
            j->jpn*n->rz*n->sy + j->jpr*n->ny*n->sz - j->jpn*n->ry*n->sz;
  */
}

/* 3D jump condition for pressure gradient
jpx =   j->jps*n->nz*n->ry - j->jps*n->ny*n->rz - j->jpr*n->nz*n->sy + j->jpn*n->rz*n->sy + j->jpr*n->ny*n->sz - j->jpn*n->ry*n->sz;
jpy =-(j->jps*n->nz*n->rx) + j->jps*n->nx*n->rz + j->jpr*n->nz*n->sx - j->jpn*n->rz*n->sx - j->jpr*n->nx*n->sz + j->jpn*n->rx*n->sz;
jpz =   j->jps*n->ny*n->rx - j->jps*n->nx*n->ry - j->jpr*n->ny*n->sx + j->jpn*n->ry*n->sx + j->jpr*n->nx*n->sy - j->jpn*n->rx*n->sy; 
det =    n->nz*n->ry*n->sx - n->ny*n->rz*n->sx  - n->nz*n->rx*n->sy  + n->nx*n->rz*n->sy  + n->ny*n->rx*n->sz  - n->nx*n->ry*n->sz;
determinant of orthongonal matrix is ONE! 
TODO: use 'det' as test to make sure it is '+1' and not '-1'
*/
