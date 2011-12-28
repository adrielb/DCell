#include "DWorld.h"

#undef __FUNCT__
#define __FUNCT__ "DWorldCreate"
PetscErrorCode DWorldCreate( Reaction rxn,
    iCoor num, Coor len, DWorld *world)
{
  struct _DWorld *w;
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  ierr = PetscNew(struct _DWorld, &w); CHKERRQ(ierr);
  if( num.z < 2 )
  {
    ierr = DACreate2d(PETSC_COMM_WORLD,
        DA_NONPERIODIC, DA_STENCIL_STAR,
        num.x,num.y,
        PETSC_DECIDE,PETSC_DECIDE,
        rxn->dof,1,0,0,&w->da); CHKERRQ(ierr);
  } else {
    ierr = DACreate3d(PETSC_COMM_WORLD,//MPI Communicator   
      DA_NONPERIODIC,   //DA_NONPERIODIC, DA_XPERIODIC, DA_YPERIODIC, DA_XYPERIODIC
      DA_STENCIL_STAR,  //DA_STENCIL_BOX or DA_STENCIL_STAR
      num.x,num.y,num.z,//Global array dimension
      PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE,//Number procs per dim
      rxn->dof,    //Number of chemical species
      1,           //stencil width
      0,0,0,       //specific array of nodes
      &w->da); CHKERRQ(ierr);
  }
  w->len = len;
  w->rxnWorld = rxn;
  ierr = DAGetMatrix(w->da, MATMPIAIJ, &w->J); CHKERRQ(ierr);
  ierr = DAGetMatrix(w->da, MATMPIAIJ, &w->L); CHKERRQ(ierr);
  MatSetFromOptions(w->J);
  
  DALocalInfo info;
  DAGetLocalInfo(w->da,&info);
  //TODO: sync dx with FluidField, so that only one 'dx' is ever defined
  PetscReal dx = len.x / (num.x-2),
            dy = len.y / (num.y-2),
            dz = len.z / (num.z-2);
  AssembleDiffusionCartesian(info, w->L,rxn->D,dx*dx,dy*dy,dz*dz);
  
  ierr = MatAssemblyBegin(w->L,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd  (w->L,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatCopy(w->L,w->J,SAME_NONZERO_PATTERN); CHKERRQ(ierr);
  
  ierr = DACreateGlobalVector(w->da, &w->global); CHKERRQ(ierr);
  
  *world = w;
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DWorldDestroy"
PetscErrorCode DWorldDestroy(DWorld world)
{
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  ierr = ReactionDestroy(world->rxnWorld); CHKERRQ(ierr);
  ierr = DADestroy(world->da); CHKERRQ(ierr);
  ierr = VecDestroy(world->global); CHKERRQ(ierr);
  ierr = MatDestroy(world->J); CHKERRQ(ierr);
  ierr = MatDestroy(world->L); CHKERRQ(ierr);
  ierr = PetscFree(world); CHKERRQ(ierr);
  
  PetscFunctionReturn(0);
}



#undef __FUNCT__
#define __FUNCT__ "Monitor"
PetscInt EVENT_Monitor;
PetscErrorCode Monitor( TS ts,PetscInt step,PetscReal time,Vec global,void *ctx)
{
  char *format = "/home/abergman/Research/DCell/temp/vec.%d.Real64";
  char name[256];
  PetscViewer binv;
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  PetscLogEventBegin(EVENT_Monitor,0,0,0,0);
//  PetscLogEventRegister(&EVENT_Monitor,"Monitor", 0);
  
//  if( step%10 != 0 ) 
//    PetscFunctionReturn(0);

  sprintf(name, format, step);
  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,name,FILE_MODE_WRITE,&binv); CHKERRQ(ierr);
  ierr = VecView(global, binv); CHKERRQ(ierr);
  ierr = PetscViewerDestroy(binv); CHKERRQ(ierr);
  
  PetscLogEventEnd(EVENT_Monitor,0,0,0,0);
  PetscFunctionReturn(0);
}
//PetscLogEventRegister(&EVENT_RHSFunction,"RHSFunction", 0);
//PetscLogEventRegister(&EVENT_RHSJacobian,"RHSJacobian", 0);