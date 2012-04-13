#include "LevelSetMethod.h"

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc, char **args) {

  PetscErrorCode ierr;
  ierr = DCellInit(); CHKERRQ(ierr);

  PetscReal dx = 1./513;
  Coor   dh = {dx,dx,0};
  Coor  len = {4,2,0};
  iCoor pos = {0,0,0};
  iCoor size = {len.x/dx,len.y/dx,0};

  LevelSet ls;
  ierr = LevelSetCreate(dh,pos,size,&ls); CHKERRQ(ierr);
  Grid g = ls->phi;
  ierr = VecSet(g->v, -1); CHKERRQ(ierr);

  //Draw border
  int b,BUF = 3;
  PetscReal **phi;
  ierr = GridGet(g,&phi); CHKERRQ(ierr);
  for (int i = 0; i < size.x; ++i) {
    for (b = 0; b < BUF; ++b) {
      phi[         b][i] = 1;
      phi[size.y-1-b][i] = 1;
    }
  }
  for (int j = 0; j < size.y; ++j) {
    for (b = 0; b < BUF; ++b) {
      phi[j][         b] = 1;
      phi[j][size.x-1-b] = 1;
    }
  }

  //Draw posts
  PetscReal w = 0.1,
            h = 0.3;
  for (PetscReal c = w; c < len.x; c+=1) {
    Coor lo = {c,  0,0};
    Coor hi = {c+w,h,0};
    ierr = DrawRect(g, lo, hi); CHKERRQ(ierr);
  }

  ierr = GridWrite(ls->phi,0); CHKERRQ(ierr);
  ierr = LevelSetInitializeFromImage(ls); CHKERRQ(ierr);
  ierr = GridWrite(ls->phi,1); CHKERRQ(ierr);
  ierr = LevelSetWriteIrregularNodeList(ls->irregularNodes,0); CHKERRQ(ierr);

  ierr = LevelSetDestroy(ls); CHKERRQ(ierr);
	ierr = DCellFinalize(); CHKERRQ(ierr);
  return 0;
}

int DrawRect( Grid g, Coor lo, Coor hi )
{
  Coor dh = g->d;
  iCoor p, q;
  int i,j,t;
  PetscErrorCode ierr;
  ierr = GridGetBounds(g,&p,&q); CHKERRQ(ierr);
  t = lo.x/dh.x;
  p.x = t < p.x ? p.x : t;
  t = lo.y/dh.y;
  p.y = t < p.y ? p.y : t;

  t = hi.x/dh.x;
  q.x = q.x < t ? q.x : t;
  t = hi.y/dh.y;
  q.y = q.y < t ? q.y : t;

  PetscReal **phi;
  ierr = GridGet(g,&phi); CHKERRQ(ierr);
  for (j = p.y; j < q.y; ++j) {
    for (i = p.x; i < q.x; ++i) {
      phi[j][i] = 1;
    }
  }
  return 0;
}
