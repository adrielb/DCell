#include "Utilities.h"
#include "petscda.h"

#undef __FUNCT__
#define __FUNCT__ "DomainSpecsCreateDefaultInterior"
PetscErrorCode DomainSpecsCreateDefaultInterior( DA da, DomainSpecs *ds )
{
  DomainSpecs spec;
  DALocalInfo info;
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  ierr = PetscNew(struct _DomainSpecs, &spec); CHKERRQ(ierr);
  ierr = DAGetLocalInfo(da, &info); CHKERRQ(ierr);
  
  spec->s.x = 2;
  spec->s.y = 2;
  spec->s.z = 2;
  
  int g  = 2;
  spec->e.x = info.mx-1-g;
  spec->e.y = info.my-1-g;
  spec->e.z = info.mz-1-g;
  
  *ds = spec;
  PetscFunctionReturn(0);
}

PetscTruth DomainSpecsIsInterior( DomainSpecs ds, iCoor g )
{
  iCoor s = ds->s, e = ds->e;
  if( s.x <= g.x && g.x <= e.x &&
      s.y <= g.y && g.y <= e.y &&
      s.z <= g.z && g.z <= g.z )
    return PETSC_TRUE;
  return PETSC_FALSE;
}

/*
PetscInt DomainSpecsToGlobalIndex( DomainSpecs, iCoor g)
{
  return 
}*/

/* Build index set of bc nodes local to this processor
 * int bc[] = (e.x-s.x)*(e.y-s.y)*(e.z-s.z);
 * for( x, y, z)
 *   if( exterior )
 *     bc[count] = (x,y,z) //global indexes
 * then shrink array
 * 
 * VecSetValues( 
 */