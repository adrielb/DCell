#include "ImmersedInterfaceMethod.h"
#include "IIM_private.h"

void InterfacialForceAdhesion( IIMIrregularNode *n, void *context );
void InterfacialForceCurvature( IIMIrregularNode *n, void *context );

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc, char **args)
{
  PetscErrorCode  ierr;
  PetscFunctionBegin;
  ierr = DCellInit(); CHKERRQ(ierr);

  PetscReal dx = 0.2;
  Coor dh = {dx, dx, dx};
  Coor center = { 2, 2, 0};
  PetscReal radius = 1;

  LevelSet ls;
  ierr = LevelSetInitializeToCircle(dh, center, radius, &ls); CHKERRQ(ierr);
//  ierr = LevelSetInitializeToStar3D(dh, center, radius, 0.4, 5, &ls); CHKERRQ(ierr);
  ierr = GridWrite(ls->phi,0); CHKERRQ(ierr);

  IIM iim;
  ierr = IIMCreate( ls->phi->is2D, &iim); CHKERRQ(ierr);
  ierr = IIMSetForceComponents( iim, InterfacialForceCurvature); CHKERRQ(ierr);
  ierr = IIMSetForceContext( iim, &dx); CHKERRQ(ierr);

  const Coor  f = (Coor){-0.1,  -0.1, 0};
  const Coor df = (Coor){ 0.06, 0.06, 1};
  ierr = IIMSetFluidCoordinateSystem( iim, f, df); CHKERRQ(ierr);

  ierr = IIMUpdateSurfaceQuantities( iim, ls); CHKERRQ(ierr);
//  ierr = SpatialIndexPrint(iim->sidx); CHKERRQ(ierr);
  ierr = IIMWriteIrregularNodeList( iim, "phi", 0); CHKERRQ(ierr);
  ierr = ArrayWrite( ls->band, 0); CHKERRQ(ierr);

  ierr = LevelSetDestroy(ls); CHKERRQ(ierr);
  ierr = IIMDestroy(iim); CHKERRQ(ierr);
  
  ierr = DCellFinalize(); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

void InterfacialForceAdhesion( IIMIrregularNode *n, void *context )
{
  const int gx = 0, gy = -1;
  n->F1 = n->nx * gx + n->ny * gy;
  n->F2 = -n->ny * gx + n->nx * gy;
}

void InterfacialForceCurvature( IIMIrregularNode *n, void *context )
{
  PetscReal dx = *(PetscReal*)context;
  n->F1 = n->X.y * dx;
  n->F2 = 0;
  n->F3 = 0;
}
