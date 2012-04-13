#include "MicrofluidicSimulator.h"

PetscErrorCode InitializeVectors( UserContext* );
PetscErrorCode AssemblePressureMatrx( UserContext* );
PetscErrorCode AssemblePressureRHS( UserContext *uc );
PetscErrorCode SolvePressure( UserContext *uc );
PetscErrorCode AssembleVelocityMatrix( UserContext* );
PetscErrorCode AssembleVelocityRHS( UserContext *uc  );
PetscErrorCode ComputeShearStress(UserContext *uc);
PetscErrorCode ConservationTest(UserContext *uc);

#undef __FUNCT__
#define __FUNCT__ "PressurePoisson"
PetscInt EVENT_PressurePoisson;
PetscErrorCode PressurePoisson( UserContext *uc )
{
  Mat A;
  Vec b, px, py, p, u, v, ss, c;
  KSP ksp;
  PetscErrorCode ierr;
  
  
  PetscFunctionBegin;
  PetscLogEventBegin(EVENT_PressurePoisson,0,0,0,0);
  PetscLogEventRegister(&EVENT_PressurePoisson,"PressurePoisson", 0);
  
  /*  Assemble and Solve for pressure */
ierr = PetscPrintf(PETSC_COMM_WORLD, "***********************\nPRESSURE\n"); CHKERRQ(ierr);
  ierr = AssemblePressureMatrx( uc ); CHKERRQ(ierr);
  ierr = AssemblePressureRHS( uc ); CHKERRQ(ierr);
  ierr = SolvePressure( uc ); CHKERRQ(ierr);
  ierr = WriteResult( uc, uc->p, "p_pressure" ); CHKERRQ(ierr);
    
  /*  Assemble and Solve for Velocity  */
ierr = PetscPrintf(PETSC_COMM_WORLD, "***********************\nVELOCITY\n"); CHKERRQ(ierr);
  ierr = AssembleVelocityRHS( uc ); CHKERRQ(ierr);
  ierr = AssembleVelocityMatrix( uc ); CHKERRQ(ierr);
  
  ierr = KSPSetOperators(uc->ksp,uc->A, uc->A, SAME_PRECONDITIONER); CHKERRQ(ierr);
  ierr = KSPSolve(uc->ksp,uc->px,uc->u);CHKERRQ(ierr);
  ierr = KSPSolve(uc->ksp,uc->py,uc->v);CHKERRQ(ierr);
  ierr = WriteResult( uc, uc->u, "u_velocity" ); CHKERRQ(ierr);
  ierr = WriteResult( uc, uc->v, "v_velocity" ); CHKERRQ(ierr);

  ierr = ComputeShearStress(uc); CHKERRQ(ierr);
  ierr = WriteResult( uc, uc->ss,"shear_stress" ); CHKERRQ(ierr);
  ierr = ConservationTest(uc); CHKERRQ(ierr);
  
  /*  Output indexing  */
  PetscReal *idx;
  VecGetArray(uc->b,&idx);
  for( int i = 0; i < uc->numNodes; i++)
    idx[i] = i;
  VecRestoreArray(uc->b,&idx);
  WriteResult(uc, uc->b, "indexes");
  VecSet(uc->b, 0.);
  
  ierr = VecDestroy(uc->c); CHKERRQ(ierr);
  ierr = VecDestroy(uc->ss); CHKERRQ(ierr);
  ierr = MatDestroy(uc->A); CHKERRQ(ierr);
  ierr = KSPDestroy(uc->ksp); CHKERRQ(ierr);
  ierr = VecDestroy(uc->b); CHKERRQ(ierr);
  ierr = VecDestroy(uc->p); CHKERRQ(ierr);
  ierr = VecDestroy(uc->px); CHKERRQ(ierr);
  ierr = VecDestroy(uc->py); CHKERRQ(ierr);
  ierr = VecDestroy(uc->u); CHKERRQ(ierr);
  ierr = VecDestroy(uc->v); CHKERRQ(ierr);
  PetscLogEventEnd(EVENT_PressurePoisson,0,0,0,0);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "InitializeVectors"
PetscErrorCode InitializeVectors( UserContext* uc)
{
  PetscErrorCode ierr;
  
  
  PetscFunctionBegin;

  ierr = VecCreate(PETSC_COMM_WORLD, &uc->b); CHKERRQ(ierr);
  ierr = VecSetSizes(uc->b, uc->numNodes, uc->numNodes); CHKERRQ(ierr);
  ierr = VecSetType(uc->b, VECSEQ); CHKERRQ(ierr);
  
  ierr = VecDuplicate(uc->b,&uc->p);CHKERRQ(ierr);
  ierr = VecDuplicate(uc->b,&uc->u);CHKERRQ(ierr);
  ierr = VecDuplicate(uc->b,&uc->v);CHKERRQ(ierr);
  ierr = VecDuplicate(uc->b,&uc->px);CHKERRQ(ierr);
  ierr = VecDuplicate(uc->b,&uc->py);CHKERRQ(ierr);
  ierr = VecDuplicate(uc->b,&uc->c);CHKERRQ(ierr);
  ierr = PetscMalloc(uc->n * sizeof(PetscReal), &uc->imageResult ); CHKERRQ(ierr);
    
  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "AssemblePressureMatrx"
PetscErrorCode AssemblePressureMatrx( UserContext* uc )
{
  PetscErrorCode  ierr;
  PetscScalar val[5];
  PetscInt i, j;
  BCNode *bcn;
  Node *n;
  
  PetscFunctionBegin;
  
  ierr = MatCreate(PETSC_COMM_WORLD, &uc->A); CHKERRQ(ierr);
  ierr = MatSetSizes(uc->A, PETSC_DECIDE, PETSC_DECIDE, uc->numNodes, uc->numNodes); CHKERRQ(ierr);
  ierr = MatSetType(uc->A, MATUMFPACK); CHKERRQ(ierr);
  ierr = MatSeqAIJSetPreallocation(uc->A, 5, PETSC_NULL); CHKERRQ(ierr);
//  ierr = MatCreateSeqSBAIJ(PETSC_COMM_WORLD,1,uc->numNodes, uc->numNodes,3, PETSC_NULL, &uc->A); CHKERRQ(ierr);
//  ierr = MatSetOption(uc->A, MAT_SYMMETRIC); CHKERRQ(ierr);
//  ierr = MatSetOption(uc->A, MAT_USE_INODES); CHKERRQ(ierr);
//  ierr = MatSetOption(uc->A, MAT_IGNORE_LOWER_TRIANGULAR); CHKERRQ(ierr);
  
  for (i = 0; i < uc->numNodes; ++i)
  {
  	n = &uc->nodes[i];
  	for (j = 0; j < n->numNei; ++j)
		{
			val[j] = negone;
		}
		val[n->numNei] = n->numNei;
		ierr = MatSetValues(uc->A, 1, &i, n->numNei+1, n->nei, val, INSERT_VALUES); CHKERRQ(ierr);
	}
  
  PetscInt idx;
  
  for( i = 0; i < uc->numBC; ++i)
  {
    idx = uc->imageToNode[uc->bcToImage[i]];
    n = &uc->nodes[idx];
    for( j = 0; j < n->numNei; ++j)
    {
      ierr = MatSetValue(uc->A, idx, n->nei[j], zero, INSERT_VALUES); CHKERRQ(ierr);
      ierr = MatSetValue(uc->A, n->nei[j], idx, zero, INSERT_VALUES); CHKERRQ(ierr);
    }
    ierr = MatSetValue(uc->A, idx, idx, one, INSERT_VALUES); CHKERRQ(ierr);
  }

  ierr = MatSetOption(uc->A, MAT_NEW_NONZERO_LOCATION_ERR); CHKERRQ(ierr);  
  ierr = MatSetOption(uc->A, MAT_NO_NEW_NONZERO_LOCATIONS); CHKERRQ(ierr);
  ierr = MatAssemblyBegin(uc->A, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(uc->A, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatStoreValues(uc->A); CHKERRQ(ierr);
  

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "AssemblePressureRHS"
PetscErrorCode AssemblePressureRHS( UserContext *uc )
{
  PetscErrorCode ierr;
  PetscInt bcidx, bcnodeidx;
  PetscReal *b;
  Node   *n;
  int i,j;
  unsigned char color;
  
  PetscFunctionBegin;
  VecGetArray(uc->b, &b);

  for (i = 0; i < uc->numBC; ++i)
  {
  	bcidx = uc->bcToImage[i];
    color = uc->filedata[bcidx];
    bcnodeidx = uc->imageToNode[bcidx];
  	n 	  = &uc->nodes[bcnodeidx];
  	b[bcnodeidx]  = uc->BCLabels[color].pressureBC;
  	
    for(j = 0; j < n->numNei; ++j)
    {
    	if( uc->filedata[uc->nodes[n->nei[j]].imageIndex] == FLUIDIC_LAYER_COLOR )
    	{
    		b[n->nei[j]] += uc->BCLabels[color].pressureBC;
    	}
    }
  }
  
  VecRestoreArray(uc->b, &b);
//  VecView( uc->b, PETSC_VIEWER_STDOUT_SELF ); 
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "SolvePressure"
PetscInt EVENT_SolvePressure;
PetscErrorCode SolvePressure(UserContext *uc)
{
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  PetscLogEventBegin(EVENT_SolvePressure,0,0,0,0);

  ierr = KSPCreate(PETSC_COMM_WORLD, &uc->ksp); CHKERRQ(ierr);
  ierr = KSPSetType(uc->ksp,KSPCG); CHKERRQ(ierr);
  ierr = KSPSetOperators(uc->ksp,uc->A,uc->A,SAME_NONZERO_PATTERN); CHKERRQ(ierr);
  ierr = KSPSetFromOptions(uc->ksp); CHKERRQ(ierr);
  ierr = KSPSetTolerances(uc->ksp, 0.000001, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT); CHKERRQ(ierr);
  
  ierr = KSPSetType( uc->ksp, KSPPREONLY); CHKERRQ(ierr);
  PC pc;
  ierr = KSPGetPC(uc->ksp, &pc); CHKERRQ(ierr);
  ierr = PCSetType(pc, PCLU); CHKERRQ(ierr);
  
  ierr = KSPSolve(uc->ksp,uc->b,uc->p);CHKERRQ(ierr);
  
  KSPConvergedReason reason;
  KSPGetConvergedReason(uc->ksp,&reason);
  PetscPrintf(PETSC_COMM_WORLD,"Pressure KSPConvergedReason: %D\n", reason);
  // if( reason < 0 )    SETERRQ(PETSC_ERR_CONV_FAILED, "Exiting: Failed to converge\n");
  
  PetscLogEventEnd(EVENT_SolvePressure,0,0,0,0);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "AssembleVelocityMatrix"
PetscErrorCode AssembleVelocityMatrix( UserContext *uc )
{
  PetscErrorCode  ierr;
  PetscScalar     val[5];
  PetscReal       z2 = PetscSqr( DELTA_Z ),
                   v = (2 * PetscSqr( DELTA_X ) + 4 * z2 ) / z2;
  int i, j;
  
  PetscFunctionBegin;
  
  ierr = MatRetrieveValues(uc->A); CHKERRQ(ierr);
  ierr = MatZeroEntries(uc->A); CHKERRQ(ierr);
  
  Node *n;
  for( i = 0; i < uc->numNodes; ++i)
	{
		n = &uc->nodes[i];
		if( n->numNei == 4 )
		{
			for (j = 0; j < 4; ++j)
			{
				if(uc->nodes[n->nei[j]].numNei == 4)
				{
					val[j] = negone;
				} else {
          val[j] = zero;
        }
			}
      val[4] = v;
			ierr = MatSetValues(uc->A,1, &i, 5, n->nei, val, INSERT_VALUES); CHKERRQ(ierr);
		} else {
			ierr = MatSetValue(uc->A, i,  i, one, INSERT_VALUES); CHKERRQ(ierr);
		}
		
	}

  ierr = MatAssemblyBegin(uc->A, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(uc->A, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

//  ierr = PetscViewerSetFormat(PETSC_VIEWER_STDOUT_SELF, PETSC_VIEWER_ASCII_DENSE); CHKERRQ(ierr);
//  ierr = MatView(uc->A, PETSC_VIEWER_STDOUT_SELF); CHKERRQ(ierr);
  
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "AssembleVelocityRHS"
PetscErrorCode AssembleVelocityRHS( UserContext *uc  )
{
  PetscErrorCode ierr;
  int i,j;
  Node *n;
  PetscScalar *p, *px, *py;
  PetscInt *nei;
  PetscReal v = DELTA_X / (two * VISCOSITY );
  
  PetscFunctionBegin;
  
  VecGetArray(uc->p , &p );
  VecGetArray(uc->px, &px);
  VecGetArray(uc->py, &py);
  
  for( i = 0; i < uc->numNodes; ++i)
  {
    n = &uc->nodes[i];
    if( n->numNei == 4 )
    {
      nei = &n->nei[0];
      px[i] = (p[nei[1]] - p[nei[2]]) * v;
      py[i] = (p[nei[0]] - p[nei[3]]) * v;
    } else {
      px[i] = zero;
      py[i] = zero;
    }
  }
  
  VecRestoreArray(uc->p , &p);
  VecRestoreArray(uc->px, &px);
  VecRestoreArray(uc->py, &py);

//printf("px:\n");  
//VecView(uc->px,   PETSC_VIEWER_STDOUT_SELF );
  
  PetscFunctionReturn(0);
}



#undef __FUNCT__
#define __FUNCT__ "ComputeShearStress"
PetscInt EVENT_ComputeShearStress;
PetscErrorCode ComputeShearStress(UserContext *uc)
{
  PetscErrorCode ierr;
  Vec u2, v2;
  PetscReal scale = VISCOSITY / DELTA_Z;
  
// mu * (u^2 + v^2)^0.5 / dz  
  PetscFunctionBegin;
  PetscLogEventBegin(EVENT_ComputeShearStress,0,0,0,0);
  
  VecDuplicate(uc->v, &uc->ss);
  VecDuplicate(uc->v, &u2);
  VecDuplicate(uc->v, &v2);
  
  VecPointwiseMult( u2, uc->u, uc->u);
  VecPointwiseMult( v2, uc->v, uc->v);
  VecWAXPY(uc->ss, one, u2, v2);
  VecSqrt(uc->ss);
  VecScale(uc->ss, scale);
  
  VecDestroy(u2);
  VecDestroy(v2);
  PetscLogEventEnd(EVENT_ComputeShearStress,0,0,0,0);
  PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "ConservationTest"
PetscInt EVENT_ConservationTest;
PetscErrorCode ConservationTest(UserContext *uc)
{
  PetscReal val[4], comp, *u, *v, *b, vx, vy,u1,u2,v1,v2;
  int i,j;
  Node *n;
  PetscInt *nei;
  PetscErrorCode ierr;
printf("***************************************************\n");
  PetscFunctionBegin;
  PetscLogEventBegin(EVENT_ConservationTest,0,0,0,0);
  PetscLogEventRegister(&EVENT_ConservationTest,"ConservationTest", 0);
  VecGetArray(uc->u, &u);
  VecGetArray(uc->v, &v);
  VecGetArray(uc->b, &b);
printf("***************************************************\n");
  for (i = 0; i < uc->numNodes; ++i)
  {
    /*
    for( j = 0; j < 4; j++ ) val[j] = 0.;
    if( u[i] < 0. ) { val[2] += u[i];}
    if( u[i] > 0. ) { val[1] -= u[i];}
    if( v[i] < 0. ) { val[3] += v[i];}
    if( v[i] > 0. ) { val[0] -= v[i];}
    for( j = 0; j < 4; j++ ) comp += val[i];*/
    n = &uc->nodes[i];
    if( n->numNei == 4 )
    {
      nei = &n->nei[0];
      vx = u[nei[1]] - u[nei[2]];
      vy = v[nei[0]] - v[nei[3]];
      u1 = u[nei[1]] - u[nei[2]];
      b[i] = vx + vy;
    } else {
      b[i] = 0.;
    }
  }
  printf("***************************************************\n");
  VecRestoreArray(uc->b, &b);
  VecRestoreArray(uc->v, &v);
  VecRestoreArray(uc->u, &u);
  printf("***************************************************\n");
  WriteResult(uc, uc->b, "conservation");
  PetscLogEventEnd(EVENT_ConservationTest,0,0,0,0);
  PetscFunctionReturn(0);
}
