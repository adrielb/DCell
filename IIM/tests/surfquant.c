#include "ImmersedInterfaceMethod.h"

void InterfacialForceAdhesion( IIMIrregularNode *n, void *context );
void InterfacialForceCurvature( IIMIrregularNode *n, void *context );

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc, char **args)
{
  PetscErrorCode  ierr;
  PetscFunctionBegin;
  ierr = DCellInit(); CHKERRQ(ierr);

  PetscReal dx = 0.01;
  Coor dh = {dx, dx, dx};
  Coor center = { 0, 0, 0};
  PetscReal radius = 1;

  LevelSet ls;
  ierr = LevelSetInitializeToSphere(dh, center, radius, &ls); CHKERRQ(ierr);
//  ierr = LevelSetInitializeToStar3D(dh, center, radius, 0.4, 5, &ls); CHKERRQ(ierr);
  ierr = GridWrite(ls->phi,0); CHKERRQ(ierr);

  PetscBool is2D = dh.z == 0.0;
  IIM iim;
  ierr = IIMCreate( is2D, &iim); CHKERRQ(ierr);
  ierr = IIMSetForceComponents(iim, InterfacialForceCurvature); CHKERRQ(ierr);
  ierr = IIMSetForceContext(iim, &dx); CHKERRQ(ierr);

  ierr = IIMUpdateSurfaceQuantities( iim, ls); CHKERRQ(ierr);
//  ierr = SpatialIndexPrint(iim->sidx); CHKERRQ(ierr);
  ierr = LevelSetWriteIrregularNodeList( ls, 0); CHKERRQ(ierr);
  ierr = ArrayWrite(ls->band,"band",0); CHKERRQ(ierr);

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
