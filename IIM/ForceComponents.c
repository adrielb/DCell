#include "ImmersedInterfaceMethod.h"
#include "IIM_private.h"

void InterfacialForceSurfaceTension( IIMIrregularNode *n, void *context )
{
  n->F1 = -n->k;
  n->F2 = 0;
}

void InterfacialForceElastic( IIMIrregularNode *n, void *context )
{
  n->F1 = 0;
  n->F2 = 0;
}

void InterfacialForceSimple( IIMIrregularNode *n, void *context )
{
  n->F1 = 10;
  n->F2 = 0;
}
