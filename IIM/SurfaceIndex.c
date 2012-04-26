#include "ImmersedInterfaceMethod.h"

#undef __FUNCT__
#define __FUNCT__ "IIMUpdateSurfaceIndex_2D"
PetscErrorCode IIMUpdateSurfaceIndex_2D( IIM iim, LevelSet ls )
{
  int i;
  int len = ArrayLength(ls->irregularNodes);
  IrregularNode *nodes = ArrayGetData(ls->irregularNodes);
  iCoor p,q;
  Coor lo, hi;
  Coor dh = { 1., 1., 1.};
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = GridGetBounds(ls->phi, &p, &q); CHKERRQ(ierr);
  lo.x = p.x;
  lo.y = p.y;
  lo.z = p.z;
  hi.x = q.x;
  hi.y = q.y;
  hi.z = q.z;
  ierr = SpatialIndexClear( iim->sidx ); CHKERRQ(ierr);
  ierr = SpatialIndexSetDomain( iim->sidx, lo, hi, dh ); CHKERRQ(ierr);
  for ( i = 0; i < len; ++i) {
    ierr = SpatialIndexInsertPoint( iim->sidx, nodes[i].X, &nodes[i] ); CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "IIMUpdateSurfaceIndex_3D"
PetscErrorCode IIMUpdateSurfaceIndex_3D( IIM iim, LevelSet ls )
{
  int i;
  int len = ArrayLength(ls->irregularNodes);
  IrregularNode *nodes = ArrayGetData(ls->irregularNodes);
  iCoor p,q;
  Coor lo, hi;
  Coor dh = { 1., 1., 1.};
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = GridGetBounds(ls->phi, &p, &q); CHKERRQ(ierr);
  lo.x = p.x;
  lo.y = p.y;
  lo.z = p.z;
  hi.x = q.x;
  hi.y = q.y;
  hi.z = q.z;
  ierr = SpatialIndexClear( iim->sidx ); CHKERRQ(ierr);
  ierr = SpatialIndexSetDomain( iim->sidx, lo, hi, dh ); CHKERRQ(ierr);
  for ( i = 0; i < len; ++i) {
    ierr = SpatialIndexInsertPoint( iim->sidx, nodes[i].X, &nodes[i] ); CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}
