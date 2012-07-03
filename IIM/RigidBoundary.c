#include "ImmersedInterfaceMethod.h"

IIMRigidBoundaryCreate( IIM iim, KSP ksp, LevelSet ls ) {
  // Setup solver
  //TODO: use superlu for parallel rigid boundary layout
  ierr = KSPCreate(PETSC_COMM_WORLD,&iim->kspRB); CHKERRQ(ierr);
  ierr = KSPSetType(iim->kspRB, KSPPREONLY); CHKERRQ(ierr);
  PC pc;
  ierr = KSPGetPC(iim->ksp, pc); CHKERRQ(ierr);
  ierr = PCSetType(pc, PCLU); CHKERRQ(ierr);

  //TODO: count the number of rigid boundary forces 'nb'

  ierr = IrregularNodeListUpdate(center,ls); CHKERRQ(ierr);
  ierr = IrregularNodeListUpdate(uface,ls); CHKERRQ(ierr);
  ierr = IrregularNodeListUpdate(vface,ls); CHKERRQ(ierr);

  int nb = 10; // Matrix A of size nb x nb
  ierr = MatCreateSeqDense(MPI_COMM_WORLD,nb,nb,PETSC_NULL,&iim->A); CHKERRQ(ierr);

  ierr = VecCreateSeq(PETSC_COMM_SELF,nb,&iim->fb); CHKERRQ(ierr);
  ierr = VecDuplicate(iim->fb,&iim->b); CHKERRQ(ierr);

  int *rows;
  ierr = PetscMalloc(sizeof(int)*nb,&rows); CHKERRQ(ierr);
  for (int i = 0; i < nb; ++i) {
    rows[i] = i;
  }

  for (int i = 0; i < nb; ++i) {
    // set fb = [0,...,0,1,0,...0], so fb[i] = 1
    ierr = PetscMemzero(rhs,sizeof(PetscReal)*rhs_len); CHKERRQ(ierr);
    rhs[ifb[i]] = 1;
    // solve the fluid eq
    ierr = KSPSolve(ksp,rhs,vel); CHKERRQ(ierr);
    // interpolate U at the rigid boundary to obtain column for A

    // set column of A
    ierr = MatSetValues(iim->A,nb,rows,1,&i,vals,INSERT_VALUES); CHKERRQ(ierr);
  }

}

IIMRigidBoundarySolve( IIM iim, LevelSet ls ) {

  // Interpolate velocity at rigid boundary
  ierr = IIMInterfaceVelocity(iim,da,vecvel,ga,ls); CHKERRQ(ierr);

  // Solve [A]fb = U0
  ierr = KSPSolve(kspRB,U0,fb); CHKERRQ(ierr);
}

void InterfacialForceSurfaceTension( IrregularNode *n )
{
  n->F1 = 1;
  n->F2 = 0;
}
