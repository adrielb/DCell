#include "ImmersedInterfaceMethod.h"
#include "IIM_private.h"

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc, char **args)
{
  PetscErrorCode  ierr;
  PetscFunctionBegin;
  ierr = DCellInit(); CHKERRQ(ierr);

  const Coor dg = {0.5,0.5,0.5};
  const iCoor size = {4,4,0};
  const iCoor pos = {0,0,0};
  LevelSet ls;
  ierr = LevelSetCreate( dg, pos, size, &ls); CHKERRQ(ierr);
  PetscReal **h;
  ierr = GridGet(ls->phi, &h); CHKERRQ(ierr);
  h[0][0] =-1; h[0][1] = -1; h[0][2] =-1; h[0][3] =-1;
  h[1][0] =-1; h[1][1] =  1; h[1][2] = 1; h[1][3] =-1;
  h[2][0] =-1; h[2][1] = -1; h[2][2] = 2; h[2][3] =-1;
  h[3][0] =-1; h[3][1] = -1; h[3][2] =-1; h[3][3] =-1;
  ierr = GridWrite(ls->phi, 0); CHKERRQ(ierr);
  ls->phi->Interpolate = GridInterpolate_Cubic;

  IIM iim;
  ierr = IIMCreate( ls->phi->is2D, &iim); CHKERRQ(ierr);

  const Coor  f = (Coor){-0.1,  -0.1, 0};
  const Coor df = (Coor){ 0.06, 0.06, 1};
  ierr = IIMSetFluidCoordinateSystem( iim, f, df); CHKERRQ(ierr);
  ierr = IIMUpdateIrregularNodes( iim, ls ); CHKERRQ(ierr);

  ierr = IIMWriteIrregularNodeList( iim, ls->phi->name, 0 ); CHKERRQ(ierr);
  const int len = ArrayLength( iim->irregularNodes );
  printf("len: %d\n", len );

  ierr = IIMDestroy( iim ); CHKERRQ(ierr);
  ierr = LevelSetDestroy( ls ); CHKERRQ(ierr);
  ierr = DCellFinalize(); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


