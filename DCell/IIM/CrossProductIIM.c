#include "ImmersedInterfaceMethod.h"

// mu L(u) + px == mu uc - pxc
//  div(u)      ==   -uc

#undef __FUNCT__
#define __FUNCT__ "IIMLaplaceCorrection"
PetscErrorCode IIMLaplaceCorrection( IIM iim, IrregularNode *n, Jump j )
{
  const PetscReal dh = ((PetscReal*)&iim->dh.x)[n->axis];
  const PetscReal dd = dh*dh;
  const PetscReal mu = *iim->mu;
  PetscReal h;
  iCoor pos;
  int *ipos = &pos.x;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  if( n->shift == CELL_CENTER ) {
    int dof = n->axis + U_FACE;
    if( n->signFace * n->signCenter > 0 ) {
      h = n->d - 1.5; // h-
    } else {
      h = 0.5 + n->d; // h+
    }
    ierr = IIMCorrection(iim, n->pos, j, n->axis, dof, n->signFace, h*dh, dd, mu ); CHKERRQ(ierr);

    if( n->signFace * n->signCenter < 0 ) {
      pos = n->pos;
      ipos[n->axis]++;
      h = n->d - 0.5; // h-
      ierr = IIMCorrection(iim, pos, j, n->axis, dof, n->signCenter, h*dh, dd, mu ); CHKERRQ(ierr);
    } else {
      pos = n->pos;
      ipos[n->axis]--;
      h = n->d - 0.5; // h+
      ierr = IIMCorrection(iim, pos, j, n->axis, dof, -n->signCenter, h*dh, dd, mu ); CHKERRQ(ierr);
    }
  } else if ( U_FACE <= n->shift && n->shift <= W_FACE){
    h = n->d - 1; // h-
    ierr = IIMCorrection(iim, n->pos, j, n->axis, n->shift, n->signCenter, h*dh, dd, mu ); CHKERRQ(ierr);
    pos = n->pos;
    ipos[n->axis]--;
    h = n->d; // h+
    ierr = IIMCorrection(iim, pos, j, n->axis, n->shift, -n->signCenter, h*dh, dd, mu ); CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "IIMVelocityGradientCorrection"
PetscErrorCode IIMVelocityGradientCorrection( IIM iim, IrregularNode *n, Jump j )
{
  PetscReal h;
  const PetscReal dh = ((PetscReal*)&iim->dh.x)[n->axis];
  iCoor pos = n->pos;
  int sign;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  if( n->shift != CELL_CENTER )
    PetscFunctionReturn(0);

  if( n->signFace * n->signCenter >= 0 ) {
    (((int*)&pos.x)[n->axis])--;
    sign = n->signFace;
  } else {
    sign = n->signCenter;
  }

  h = n->d - 0.5;
  h *= dh;

  ierr = IIMCorrection( iim, pos, j, n->axis, CELL_CENTER, sign, h, dh, 1 ); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "IIMPressureGradientCorrection"
PetscErrorCode IIMPressureGradientCorrection( IIM iim, IrregularNode *n, Jump j )
{
  PetscReal h;
  const PetscReal dh = ((PetscReal*)&iim->dh.x)[n->axis];
  PetscErrorCode ierr;

  PetscFunctionBegin;
  if( n->shift != CELL_CENTER )
    PetscFunctionReturn(0);

  // h- = x_i - a and h+ = x_i+1 - a
  h = n->d > 0.5 ? n->d - 1 : n->d;
  h *= dh;

  ierr = IIMCorrection(iim, n->pos, j, n->axis, n->axis+U_FACE, n->signCenter, h, dh, 1); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "IIMCorrection"
PetscErrorCode IIMCorrection( IIM iim, iCoor x, Jump j, int axis, VelFace dof, int sign, PetscReal h, PetscReal dd, PetscReal mu )
{
  int *coor;
  PetscReal *C;
  PetscReal *jpq  = &j.x;
  PetscReal *jpqq = &j.xx;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = ArrayAppend(iim->val, &C); CHKERRQ(ierr);
  *C = mu * sign * ( j.j + h * jpq[axis] + 0.5 * h * h * jpqq[axis] ) / dd;

  ierr = ArrayAppend(iim->coor,&coor); CHKERRQ(ierr);
  coor[0] = dof;
  coor[1] = x.x;
  coor[2] = x.y;
  coor[3] = x.z;

  PetscFunctionReturn(0);
}

void JumpPressure( IrregularNode *n, Jump *j )
{
  j->j = n->f1; //#
  j->e = n->f2_n + n->f3_t; //#
  j->n = n->f1_n; //#
  j->t = n->f1_t;
  j->nn = n->f1_nn - n->k_nn * j->e; //#
  j->tt = n->f1_tt - n->k_tt * j->e;
  j->nt = n->f1_nt - n->k_nt * j->e;
  j->ne = n->f2_nn + n->f3_nt + j->n * n->k_nn + j->t * n->k_nt; //#
  j->te = n->f2_nt + n->f3_tt + j->n * n->k_nt + j->t * n->k_tt;
  j->ee = -j->nn - j->tt; //#

  IIMLocalToGlobal_1st( n, j );
  IIMLocalToGlobal_2nd( n, j );
}

void JumpVelocity( PetscReal mu, IrregularNode *n, Jump *j, int i )
{
//  PetscReal *ft = &n->ftx;
  PetscReal *e = &n->nx, *s = &n->sx, *r = &n->rx;
// #: verified
  j->j = 0; //#
  j->n = 0; //#
  j->t = 0;
  j->e = -n->f2*s[i] / mu; //ft[i] / mu; // Ft = F - Fn.n //#
  j->nn = -n->k_nn * j->e; //#
  j->tt = -n->k_tt * j->e / mu;
  j->nt = -n->k_nt * j->e / mu;
  j->te = 0; //(j->e)_t / mu; // n->ftx_t
  j->ne = (-n->f2_n*s[i] - n->k_nn*n->f2*e[i])/mu; //# //(j->e)_n / mu; // TODO: tangential force derivatives
  PetscReal pe = n->f2_n, // + n->f3_t,
            pn = n->f1_n,
            pt = n->f1_t;
  j->ee = -j->nn - j->tt + ( pe * e[i] + pn * s[i] + pt * r[i] ) / mu; //#

  IIMLocalToGlobal_1st( n, j );
  IIMLocalToGlobal_2nd( n, j );
}
