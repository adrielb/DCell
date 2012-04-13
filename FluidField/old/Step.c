#include "FluidField.h"

/*
 * 1) Solve Laplace(p) = pc + 4/3 mu/rho Laplace(s)
 * 2) Solve Laplace(u) = px/mu + pxc - uc - 1/3rho sx
 * 3) Solve Laplace(v) = py/mu + pyc - vc - 1/3rho sy
 * *) Solve Laplace(p) = pc + div (do i need this projection step?)
 */

#undef __FUNCT__
#define __FUNCT__ "FluidFieldStep"
PetscLogEvent EVENT_FluidFieldStep;
PetscErrorCode FluidFieldStep(FluidField f, IIM iim, int lsLEN, LevelSet *ls)
{
  const iCoor X_FACE = {-1, 0, 0}, 
              Y_FACE = { 0,-1, 0},
              Z_FACE = { 0 ,0,-1},
              FACES[3] = {X_FACE,Y_FACE,Z_FACE},
              CELL_CENTER = {0,0,0};
  Mat *mat = &f->matU;
  Vec *dp  = &f->px;
  Vec *vel = &f->u;
  char *dpname[3] = {"px", "py", "pz" };
  char *dpcname[3]= {"pxc","pyc","pzc"};
  char *velname[3]= {"u",  "v",  "w"  };
  PetscReal sink;
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  PetscLogEventBegin(EVENT_FluidFieldStep,0,0,0,0);
//  PetscLogEventRegister("FluidFieldStep", COOKIE_FluidField, &EVENT_FluidFieldStep);

  //update source vector
  ierr = VecAssemblyBegin(f->source); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(f->source); CHKERRQ(ierr);
  ierr = VecSum(f->source,&sink); CHKERRQ(ierr); // Ensure mass conservation: sources + sinks = 0
  int rank;
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
  if( rank == 0 )
  {
    iCoor sinkCoor = {4,4,0};
    int idx;
    PetscReal *s;
    ierr = VecGetArray(f->source, &s); CHKERRQ(ierr);
    int a;
    for (a = 0; a < 20; ++a) {
      sinkCoor.x = a + 4;
      idx = DACoorToGlobalIndex( f->da, sinkCoor );
      s[idx] -= sink / 20.;
    }
    ierr = VecRestoreArray(f->source, &s); CHKERRQ(ierr);
  }
//  ierr = MatMult(f->matP,f->source,f->div); CHKERRQ(ierr); // 'div' used as pressure poisson rhs
  ierr = VecCopy(f->source,f->div); CHKERRQ(ierr);
  ierr = VecScale(f->div,4*f->mu/(3*f->rho)); CHKERRQ(ierr);
  
/*
 *  Laplace(p) = pc + 4/3 mu/rho Laplace(s)
 * 
 */
  ierr = VecZeroEntries(f->div); CHKERRQ(ierr);
  int i;
  for( i = 0; i < lsLEN; ++i)
  {
    ierr = IIMUpdateSurfaceQuantities(iim,CELL_CENTER,ls[i]); CHKERRQ(ierr);
    ierr = IIMLaplaceCorrection(iim, ls[i], 0, JumpPressure, f->d, f->div); CHKERRQ(ierr);
  }
  
  ierr = VecAssemblyBegin(f->div); CHKERRQ(ierr); //TODO: manually time how long this assembly takes
  ierr = VecAssemblyEnd(f->div); CHKERRQ(ierr);
  ierr = KSPSolve(f->kspP,f->div,f->p); CHKERRQ(ierr);
  ierr = VecWrite("pc", f->div); CHKERRQ(ierr);
  ierr = VecWrite("p" , f->p); CHKERRQ(ierr);

/*
 * [ px, py ] = grad(p)/mu - 1/3rho [sx, sy] + [pxc, pyc]
 * 
 */

  ierr = FluidFieldGradient( f ); CHKERRQ(ierr); // Given p & s, update px, py, pz
  int v;
  for( v = 0; v < (f->is3D ? 3 : 2); v++ )
  {
    ierr = VecWrite(dpname[v],dp[v]); CHKERRQ(ierr);
    int i;
    for( i = 0; i < lsLEN; ++i)
    {
      //TODO: combined all faces into one irreg list, so that irreg nodes = [cell center, x face, ...]
      ierr = IIMGradientCorrection(iim, ls[i], v, JumpPressure, f->d, dp[v]); CHKERRQ(ierr);
    }
  }
  for( v = 0; v < (f->is3D ? 3 : 2); v++ )
  {
    int i;
    for (i = 0; i < lsLEN; ++i)
    {
      ierr = IIMUpdateSurfaceQuantities(iim,FACES[v],ls[i]); CHKERRQ(ierr);
      ierr = IIMLaplaceCorrection( iim, ls[i], v, JumpVelocity, f->d, dp[v]); CHKERRQ(ierr);
    }
    ierr = VecAssemblyBegin(dp[v]); CHKERRQ(ierr);
  }
  
/*
 *  Laplace(u) = ( px + pxc - uc - mu/3rho sx ) / mu
 * 
 */
  for( i = 0; i < (f->is3D ? 3 : 2); i++ )
  {
    ierr = VecAssemblyEnd(dp[i]); CHKERRQ(ierr);
    ierr = VecWrite(dpcname[i],dp[i]); CHKERRQ(ierr);
    ierr = KSPSetOperators(f->kspU,mat[i],mat[i],SAME_PRECONDITIONER); CHKERRQ(ierr);
    ierr = KSPSolve(f->kspU,dp[i],vel[i]); CHKERRQ(ierr);
    ierr = VecWrite(velname[i], vel[i]); CHKERRQ(ierr);
  }
  
  
// div = du/dx + dv/dy
  ierr = FluidFieldDivergence( f ); CHKERRQ(ierr); // Update f->div from u,v,w
VecWrite("div",f->div);
  

   /* 
//Solve Laplace(phi) = div
  ierr = KSPSolve(f->kspPhi,f->div,f->p); CHKERRQ(ierr);
  ierr = FluidFieldGradient( f ); CHKERRQ(ierr); 
  
//Velocity Correction
// [u, v] -= [px, py]
  ierr = VecAXPY(f->u,-1,f->px); CHKERRQ(ierr);
  ierr = VecAXPY(f->v,-1,f->py); CHKERRQ(ierr);

  ierr = FluidFieldDivergence( f ); CHKERRQ(ierr);
  VecWrite("divZ",f->div); 
*/


  PetscLogEventEnd(EVENT_FluidFieldStep,0,0,0,0);
  PetscFunctionReturn(0);
}