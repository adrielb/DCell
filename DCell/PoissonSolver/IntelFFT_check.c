#include "IntelFFT.h"
#include "FluidField.h"

int main(int argc, char **args)
{
  PC pc;
  FluidField ff;
  DALocalInfo g;
  IS *is;
  PetscErrorCode  ierr;
  ierr = PetscInitialize(&argc, &args, (char *) 0, ""); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "Starting: %s\n", __FILE__); CHKERRQ(ierr);
  
  ierr = FluidFieldCreate(&ff); CHKERRQ(ierr);
  ierr = KSPGetPC(ff->kspP, &pc); CHKERRQ(ierr);
  ierr = PCSetType(pc, PCASM); CHKERRQ(ierr);
  KSP *subksp;
  ierr = PCASMGetSubKSP(pc, 0, 0, &subksp); CHKERRQ(ierr);
  ierr = PCASMGetLocalSubdomains(pc,0,&is); CHKERRQ(ierr);
  ierr = KSPSetType(subksp[0],KSPPREONLY); CHKERRQ(ierr);
  ierr = KSPGetPC(subksp[0],&pc); CHKERRQ(ierr);
  ierr = DAGetLocalInfo(ff->da, &g); CHKERRQ(ierr);
  iCoor s = {2,2,2};
  iCoor e = {g.mx-2, g.my-2, g.mz-2};
  ierr = PCIntelFFTSetPC( pc, g, is[0], s, e ); CHKERRQ(ierr);
  
  ierr = PetscPrintf(PETSC_COMM_WORLD, "End\n"); CHKERRQ(ierr);
  ierr = PetscFinalize(); CHKERRQ(ierr);
}
