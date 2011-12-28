#include "DWorld.h"

int main(int argc, char **args)
{
  PetscErrorCode  ierr;
  ierr = PetscInitialize(&argc, &args, (char *) 0, ""); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "Start\n"); CHKERRQ(ierr);
  
  Reaction rxn;
  ierr = ReactionCreate(1,&rxn); CHKERRQ(ierr);
  rxn->D[0] = 1e-3;
  
  int n = 100;
  iCoor num = {n,n,0};
  Coor len = {1, 1, 0};
  PetscReal dx = 1./(n-2);
  Coor dh = {dx,dx,0};
  
  DWorld w;
  ierr = DWorldCreate(rxn,num,len,&w); CHKERRQ(ierr);
  ierr = MatZeroEntries(w->L); CHKERRQ(ierr);

  Vec u, v;
  ierr = VecDuplicate(w->global,&u); CHKERRQ(ierr);
  ierr = VecDuplicate(w->global,&v); CHKERRQ(ierr);
  int i,j;
  double xU,xV,yU,yV, **us,**vs, pi = 3.141592653589793;
  ierr = DAVecGetArray(w->da,u,&us); CHKERRQ(ierr);
  ierr = DAVecGetArray(w->da,v,&vs); CHKERRQ(ierr);
  DALocalInfo info;
  DAGetLocalInfo(w->da,&info);
  PetscReal xe, ye;
  xe = info.xs+info.xm;
  ye = info.ys+info.ym;
  xe = xe == info.mx ? xe - 1 : xe;
  ye = ye == info.my ? ye - 1 : ye;
  for( j = info.ys; j < ye; ++j)
  {
    for( i = info.xs; i < xe; ++i)
    {
      xU = (i-0.5) * dx;
      yU = j * dx;
      xV = i * dx;
      yV = (j-0.5) * dx;
      us[j][i] =  2*cos(pi*yU)*sin(pi*xU)*sin(pi*xU)*sin(pi*yU);
      vs[j][i] = -2*cos(pi*xV)*sin(pi*xV)*sin(pi*yV)*sin(pi*yV);
      /*us[j][i] = 0;
      vs[j][i] = 0;*/
    }
  }
  ierr = DAVecRestoreArray(w->da,u,&us); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(w->da,v,&vs); CHKERRQ(ierr);
  
  Vec ul, vl;
  ierr = DAGetLocalVector(w->da,&ul); CHKERRQ(ierr);
  ierr = DAGetLocalVector(w->da,&vl); CHKERRQ(ierr);
  ierr = DAGlobalToLocalBegin(w->da,u,INSERT_VALUES,ul); CHKERRQ(ierr);
  ierr = DAGlobalToLocalEnd  (w->da,u,INSERT_VALUES,ul); CHKERRQ(ierr);
  ierr = DAGlobalToLocalBegin(w->da,v,INSERT_VALUES,vl); CHKERRQ(ierr);
  ierr = DAGlobalToLocalEnd  (w->da,v,INSERT_VALUES,vl); CHKERRQ(ierr);
  ierr = DAVecGetArray(w->da,ul,&us); CHKERRQ(ierr);
  ierr = DAVecGetArray(w->da,vl,&vs); CHKERRQ(ierr);
  
  AssembleDiffusionCartesian(info, w->L,rxn->D,dx*dx,dx*dx,dx*dx);
  ierr = AssembleAdvectiveTransport(info,w->L,dh,us,vs); CHKERRQ(ierr);
  ierr = MatAssemblyBegin(w->L,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(w->L,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  
  Vec one;
  VecDuplicate(w->global,&one);
  VecSet(one,1);
  
  PetscReal dt = -.1*dx;
  MatScale(w->L,dt);
  MatDiagonalSet(w->L,one,ADD_VALUES);
/*
  PetscViewer view;
  PetscViewerASCIIOpen(PETSC_COMM_SELF, "/home/abergman/Research/DCell/temp/mat.dat", &view);
  MatView(w->L, view);
*/  
  KSP ksp;
  KSPCreate(PETSC_COMM_WORLD,&ksp);
  KSPSetOperators(ksp,w->L,w->L,SAME_PRECONDITIONER);
  KSPSetInitialGuessNonzero(ksp,PETSC_TRUE);
  PC pc;
  KSPGetPC(ksp, &pc);
  PCSetType(pc,PCJACOBI);
  
  WriteVec wv;
  WriteVecCreate(v,"v",&wv);
  VecShift(v,2);
  Monitor(0,0,0,v,0);
  for( int i = 0; i < 400; i++)
  {
    for(int j = 0; j < 10; j++ )
      KSPSolve(ksp,v,v);
    Monitor(0,i,0,v,0);
    printf("%d\n",i);
  }
  
  ierr = PetscPrintf(PETSC_COMM_WORLD, "End\n"); CHKERRQ(ierr);
  ierr = PetscFinalize(); CHKERRQ(ierr);
  return 0;
}