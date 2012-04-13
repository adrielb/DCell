#include "MicrofluidicSimulator.h"

PetscErrorCode AssembleCoupledStokesEq( UserContext *uc );

#undef __FUNCT__
#define __FUNCT__ "AssembleCoupledStokesEq"
PetscInt EVENT_AssembleCoupledStokesEq;
PetscErrorCode AssembleCoupledStokesEq( UserContext *uc )
{
  PetscInt nn = uc->numNodes, size = 3*nn, row, j, idx[5];
  PetscInt colu[7], colv[7], colc[5];
  Node *n;
//  const PetscReal uh2 = VISCOSITY / ( DELTA_X * DELTA_X ),
  const PetscReal uh2 = 1.5,
              valu[7] = {-uh2,-uh2,-uh2,-uh2,4*uh2,-1,1}, 
              valv[7] = {-uh2,-uh2,-uh2,-uh2,4*uh2,-1,1},
              valc[4] = {-1,-1,1,1};
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  PetscLogEventBegin(EVENT_AssembleCoupledStokesEq,0,0,0,0);
  PetscLogEventRegister(&EVENT_AssembleCoupledStokesEq,"AssembleCoupledStokesEq", 0);
  
  Mat mat;
  MatCreate(PETSC_COMM_SELF, &mat);
  MatSetSizes(mat,size,size,size,size);
  MatSetType(mat, MATSEQAIJ);
  MatSeqAIJSetPreallocation(mat, 7, PETSC_NULL);
  for( int i = 0; i < nn; i++)
  {
    n = &uc->nodes[i];
    for( j = 0; j < 5; j++)
    {
      if( n->star[j] == -1 ) idx[j] = 0; else idx[j] = 1;
      colu[j] = n->star[j];             // div grad u columns
      colv[j] = n->star[j] + nn*idx[j]; // div grad v columns
      colc[j] = n->star[j];             // ux + vy
    }
    colu[5] = n->star[3] + 2*nn*idx[3]; // -py = p0 - p3 
    colu[6] = n->star[0] + 2*nn*idx[0]; 
    colv[5] = n->star[2] + 2*nn*idx[2]; // -px = p1 - p2
    colv[6] = n->star[1] + 2*nn*idx[1];
    colc[0] += nn*idx[0];               // vy1
    colc[3] += nn*idx[3];               // vy2
    if( n->numNei != 4 )
    {
      MatSetValue(mat,i   ,i   ,1.,INSERT_VALUES); // u = 0 no-slip BC
      MatSetValue(mat,i+nn,i+nn,1.,INSERT_VALUES); // v = 0 no-slip BC
      row = i + 2*nn; // pressure index
      MatSetValue(mat,row,row,1.,INSERT_VALUES);
      if( n->numNei <= 2 ) // pressure equal to any of its neighbors
      {
        for( j = 0; j < 4; j++)
        {
          if( n->star[j] == -1 ) continue;
          MatSetValue(mat,row,2*nn+n->star[j],-1.,INSERT_VALUES);
          break;
        }
      }
      if( n->numNei == 3 ) // pressure equal to tangent node
      {
        for( j = 0; j < 4; j++)
        {
          if( n->star[j] != -1 ) continue;
printf("%d, %d\n", i, 2*nn+n->star[3-j]);
          MatSetValue(mat,row,2*nn+n->star[3-j],-1.,INSERT_VALUES);
          break;
        }
      }
    } else {
      row = i;
      MatSetValues(mat,1,&row,7,(PetscInt*)&colu,(PetscReal*)&valu,INSERT_VALUES);
      row = i + nn;
      MatSetValues(mat,1,&row,7,(PetscInt*)&colv,(PetscReal*)&valv,INSERT_VALUES);
      row = i + 2*nn;
      MatSetValues(mat,1,&row,4,(PetscInt*)&colc,(PetscReal*)&valc,INSERT_VALUES);
    }
  }
  
  MatAssemblyBegin(mat,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(mat,MAT_FINAL_ASSEMBLY);
  
  Vec vec, res;
  VecCreateSeq(PETSC_COMM_SELF,size,&vec);
  VecDuplicate(vec,&res);
  PetscReal *v;
  VecGetArray(vec,&v);
  unsigned char color;
  PetscInt bcidx, bcnodeidx;
  for( int i = 0; i < uc->numBC; ++i)
  {
    bcidx = uc->bcToImage[i];
    color = uc->filedata[bcidx];
    bcnodeidx = uc->imageToNode[bcidx] + 2 * nn;
    printf("bcnodeidx: %d, %d, %f\n", bcnodeidx, color,uc->BCLabels[color].pressureBC);
    v[bcnodeidx] = uc->BCLabels[color].pressureBC;
    MatZeroRows(mat,1,&bcnodeidx,1.);
  }
  VecRestoreArray(vec,&v);

MatView(mat,PETSC_VIEWER_STDOUT_SELF);
  
  PetscReal *r, *b;
  VecGetArray(res,&r);
  VecGetArray(uc->b,&b);
  VecGetArray(vec,&v);
  PetscMemcpy(b, &v[2*nn], nn*sizeof(PetscReal) );
  ierr = WriteResult(uc, uc->b, "p_ic"); CHKERRQ(ierr);
  
  KSPCreate(PETSC_COMM_SELF,&uc->ksp);
  ierr = KSPSetType( uc->ksp, KSPPREONLY); CHKERRQ(ierr);
  PC pc;
  ierr = KSPGetPC(uc->ksp, &pc); CHKERRQ(ierr);
  ierr = PCSetType(pc, PCLU); CHKERRQ(ierr);
  KSPSetOperators(uc->ksp, mat, mat, DIFFERENT_NONZERO_PATTERN);
  ierr = KSPSolve(uc->ksp, vec, res); CHKERRQ(ierr);
  
  
  PetscMemcpy(b, res, nn*sizeof(PetscReal) );
  ierr = WriteResult(uc, uc->b, "u_velocity"); CHKERRQ(ierr);  
  PetscMemcpy(b, &res[nn], nn*sizeof(PetscReal) );
  ierr = WriteResult(uc, uc->b, "v_velocity"); CHKERRQ(ierr);
  PetscMemcpy(b, &res[2*nn], nn*sizeof(PetscReal) );
  ierr = WriteResult(uc, uc->b, "p_pressure"); CHKERRQ(ierr);
  
  
//  MatSetFromOptions(mat);
  
  PetscLogEventEnd(EVENT_AssembleCoupledStokesEq,0,0,0,0);
  PetscFunctionReturn(0);
}

void reg()
{
  
}