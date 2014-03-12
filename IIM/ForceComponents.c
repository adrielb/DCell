#include "ImmersedInterfaceMethod.h"

void InterfacialForceSurfaceTension( IrregularNode *n, void *context )
{
  UNUSED(context);
  n->F1 = -n->k;
  n->F2 = 0;
}

void InterfacialForceElastic( IrregularNode *n, void *context )
{
  UNUSED(context);
  n->F1 = 0;
  n->F2 = 0;
}

void InterfacialForceSimple( IrregularNode *n, void *context )
{
  UNUSED(context);
  n->F1 = 10;
  n->F2 = 0;
}
