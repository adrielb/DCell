#include "ImmersedInterfaceMethod.h"

void InterfacialForceSurfaceTension( IrregularNode *n, void *context )
{
  n->f1 = -n->k;
  n->f2 = 0;
}

void InterfacialForceElastic( IrregularNode *n, void *context )
{
  n->f1 = 0;
  n->f2 = 0;
}

void InterfacialForceSimple( IrregularNode *n, void *context )
{
  n->f1 = 10;
  n->f2 = 0;
}
