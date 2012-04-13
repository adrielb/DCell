#include "petsc.h"
#include "petscda.h"

int main( int argc, char **args)
{
  PetscErrorCode  ierr;
  ierr = PetscInitialize(&argc, &args, (char *) 0, ""); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "Start\n"); CHKERRQ(ierr);

  DA da;
  ierr = DACreate2d(PETSC_COMM_WORLD,DA_NONPERIODIC,DA_STENCIL_STAR,
		    64,64,PETSC_DECIDE,PETSC_DECIDE,1,1,0,0, &da); CHKERRQ(ierr);
  int rank;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  printf("[%d] rank\n", rank);
  Vec p;
  ierr = DACreateGlobalVector(da,&p); CHKERRQ(ierr);
  ierr = VecSet( p, 1.);
  char file[256] = "vec";
  sprintf( file, "vec-%d", rank);
  PetscViewer view;
  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,file,FILE_MODE_WRITE,&view
			       ); CHKERRQ(ierr);
  ierr = VecView(p, view); CHKERRQ(ierr);
  ierr = PetscViewerDestroy(view); CHKERRQ(ierr);  
  ierr = VecDestroy( p );
  ierr = DADestroy( da );
  ierr = PetscPrintf(PETSC_COMM_WORLD, "End\n"); CHKERRQ(ierr);
  ierr = PetscFinalize(); CHKERRQ(ierr);
}
