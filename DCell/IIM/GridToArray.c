#include "ImmersedInterfaceMethod.h"
#include "IIM_private.h"

#undef __FUNCT__
#define __FUNCT__ "IIMUpdateIrregularNodeGrid_2D"
PetscErrorCode IIMUpdateIrregularNodeGrid_2D( IIM iim, LevelSet ls )
{
  int i;
  int len = ArrayLength(ls->irregularNodes);
  IrregularNode *nodes = ArrayGetData(ls->irregularNodes);
  GridPoint *gridpoint;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscLogEventBegin(EVENT_IIMUpdateIrregularNodeGrid,0,0,0,0); CHKERRQ(ierr);
  iim->p = ls->phi->p;
  iim->n = ls->phi->n;
  ierr = ArraySetSize(iim->irregularNodeGrid,iim->n.x*iim->n.y); CHKERRQ(ierr);
  ierr = ArrayZero(iim->irregularNodeGrid); CHKERRQ(ierr);
  for ( i = 0; i < len; ++i) {
    ierr = IIMGetIrregularNodeGridPoint_2D( iim, nodes[i].pos.x, nodes[i].pos.y, &gridpoint); CHKERRQ(ierr);
    gridpoint->len++;
    if( gridpoint->idx == 0 )
      gridpoint->idx = i;
  }
  ierr = PetscLogEventEnd(EVENT_IIMUpdateIrregularNodeGrid,0,0,0,0); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "IIMUpdateIrregularNodeGrid_3D"
PetscErrorCode IIMUpdateIrregularNodeGrid_3D( IIM iim, LevelSet ls )
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  SETERRQ(PETSC_COMM_SELF,0,"NOT IMPELEMENTED");
  ierr=0;
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "IIMGetIrregularNodeGridPoint_2D"
PetscErrorCode IIMGetIrregularNodeGridPoint_2D( IIM iim, int x, int y, GridPoint **pt )
{
  int idx;
  iCoor c;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  c.x = x - iim->p.x;
  c.y = y - iim->p.y;
  idx = c.x + iim->n.x*c.y;
  ierr = ArrayGet(iim->irregularNodeGrid,idx,pt); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
