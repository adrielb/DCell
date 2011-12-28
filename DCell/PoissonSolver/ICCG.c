#include "PoissonSolver.h"
#include "petscksp.h"
#include "petscda.h"
//#include "LevelSet.h"

#undef __FUNCT__
#define __FUNCT__ "Generate2DLaplacianPeriodicBC"
PetscInt EVENT_Generate2DLaplacianPeriodicBC;
PetscErrorCode Generate2DLaplacianPeriodicBC(  PetscInt d1, PetscInt d2, Mat *mat )
{
  DA da;
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  PetscLogEventBegin(EVENT_Generate2DLaplacianPeriodicBC,0,0,0,0);
  
  ierr = DACreate2d(PETSC_COMM_SELF,//MPI Communicator   
    DA_NONPERIODIC,   //DA_NONPERIODIC, DA_XPERIODIC, DA_YPERIODIC, DA_XYPERIODIC
    DA_STENCIL_STAR, //DA_STENCIL_BOX or DA_STENCIL_STAR
    d1, d2, //Global dimension
    1, 1,          //Number procs per dim
    1,            //dof
    1,           //stencil width
    0,0,        //specific array of nodes
    &da); 
  
  
//  ierr = DAGetMatrix( da, MATSEQSBAIJ, mat); CHKERRQ(ierr);
//  ierr = MatSetOption(*mat,MAT_IGNORE_LOWER_TRIANGULAR); CHKERRQ(ierr);
  ierr = DAGetMatrix( da, MATSEQAIJ, mat); CHKERRQ(ierr);
  MatStencil row, col[5];
  PetscReal v[5] = {4, -1, -1, -1, -1}, o = 1.;
  PetscInt i,j;
  for( j = 0; j < d2; ++j)
  {
    for( i = 0; i < d1; ++i)
    {
      row.i = i;
      row.j = j;
      col[0].i = i;              col[0].j = j;
      col[1].i = i;              col[1].j = (j+1)%d2;
      col[2].i = (i+1)%d1;       col[2].j = j;
      col[3].i = (i-1+d1)%d1;    col[3].j = j;
      col[4].i = i;              col[4].j = (j-1+d2)%d2;
      ierr = MatSetValuesStencil(*mat, 1, &row, 5, col, v, INSERT_VALUES); CHKERRQ(ierr);
    }
  }
  
  i=d1/2;
  j=d2/2;
  row.i = i;
  row.j = j;
  col[0].i = i;              col[0].j = j;
  col[1].i = i;              col[1].j = (j+1)%d2;
  col[2].i = (i+1)%d1;       col[2].j = j;
  col[3].i = (i-1+d1)%d1;    col[3].j = j;
  col[4].i = i;              col[4].j = (j-1+d2)%d2;
  PetscReal vv[5] = {1, 0, 0, 0, 0};
  ierr = MatSetValuesStencil(*mat, 1, &row, 5, col, vv, INSERT_VALUES); CHKERRQ(ierr);
  
  PetscReal z = 0;
  col[0].i = i;              col[0].j = j;
  row.i    = (i-1)%d1;          row.j = j;
  ierr = MatSetValuesStencil(*mat, 1, &row, 1, col, &z, INSERT_VALUES); CHKERRQ(ierr);
  row.i    = j;                 row.j = (j-1)%d2;
  ierr = MatSetValuesStencil(*mat, 1, &row, 1, col, &z, INSERT_VALUES); CHKERRQ(ierr);
  row.i    = (i+1)%d1;          row.j = j;
  ierr = MatSetValuesStencil(*mat, 1, &row, 1, col, &z, INSERT_VALUES); CHKERRQ(ierr);
  row.i    = j;                 row.j = (j+1)%d2;
  ierr = MatSetValuesStencil(*mat, 1, &row, 1, col, &z, INSERT_VALUES); CHKERRQ(ierr);
  
  MatAssemblyBegin(*mat,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(*mat, MAT_FINAL_ASSEMBLY);
  
//  PetscInt z = 0;
//  MatZeroRows(*mat,1,&z,1.);
  
  PetscLogEventEnd(EVENT_Generate2DLaplacianPeriodicBC,0,0,0,0);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "Generate2DLapacian"
PetscInt EVENT_Generate2DLapacian;
PetscErrorCode Generate2DLapacian( PetscInt d1, PetscInt d2, Mat *m )
{
  PetscErrorCode ierr;
  DA da;
  
  PetscFunctionBegin;
  PetscLogEventBegin(EVENT_Generate2DLapacian,0,0,0,0);
  
  /* COMM, block size, num rows, num cols, num nonzeros, array nz, mat */
//  MatCreateSeqSBAIJ(PETSC_COMM_SELF, 1, ls->g2d->len, ls->g2d->len, 3, PETSC_NULL, &mat );
//  MatSetValues(Mat mat,PetscInt m,const PetscInt idxm[],PetscInt n,const PetscInt idxn[],const PetscScalar v[],InsertMode addv)
  
  ierr = DACreate2d(PETSC_COMM_SELF,//MPI Communicator   
  DA_NONPERIODIC,   //DA_NONPERIODIC, DA_XPERIODIC, DA_YPERIODIC, DA_XYPERIODIC
  DA_STENCIL_STAR, //DA_STENCIL_BOX or DA_STENCIL_STAR
  d1, d2, //Global dimension
  1, 1,          //Number procs per dim
  1,            //dof
  1,           //stencil width
  0,0,        //specific array of nodes
  &da);
  
  Mat mat;
  
  DAGetMatrix( da, MATSEQSBAIJ, &mat);

  MatStencil row, col[5];
  PetscReal v[3] = {4, -1, -1}, o = 1.;
  PetscInt i,j;
  
//  DAGetCorners(da, &xs,&ys,0, &d1,&d2,0 );
  for( j = 0; j < d2; ++j)
  {
    for( i = 0; i < d1; ++i)
    {
      row.i = i;
      row.j = j;
      if( i == 0 || j == 0 || i == d1 - 1 || j == d2 -1 )
      {
        MatSetValuesStencil(mat,1,&row,1,&row,&o, INSERT_VALUES);
      } else {
//        col[0].i = i; col[0].j = j-1;
//        col[1].i = i; col[1].j = j+1;
//        col[2].i = i; col[2].j = j;
//        col[3].i = i-1; col[3].j = j;
//        col[4].i = i+1; col[4].j = j;
        col[0].i = i; col[0].j = j;
        col[1].i = i; col[1].j = j+1;
        col[2].i = i+1; col[2].j = j;
        MatSetValuesStencil(mat, 1, &row, 3, col, v, INSERT_VALUES);
      }
    }
  }
  
  PetscReal z=0.;
  for( i = 1; i < d1-1; ++i)
  {
    row.i = i;
    row.j = d2-2;
    col[0].i = i;
    col[0].j = d2-1;
    MatSetValuesStencil(mat, 1, &row,1,col,&z, INSERT_VALUES );
  }
  for( j = 1; j < d2-1; ++j)
  {
    row.i = d1-2;
    row.j = j;
    col[0].i = d1-1;
    col[0].j = j;
    MatSetValuesStencil(mat, 1, &row,1,col,&z, INSERT_VALUES );
  }
  
  MatAssemblyBegin(mat,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(mat,MAT_FINAL_ASSEMBLY);
  *m = mat;
  
  PetscLogEventEnd(EVENT_Generate2DLapacian,0,0,0,0);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "GenerateLaplacian2DNoBC"
PetscInt EVENT_GenerateLaplacian2DNoBC;
PetscErrorCode GenerateLaplacian2DNoBC( PetscInt d1, PetscInt d2, Mat *m )
{
  PetscErrorCode ierr; 
  DA da;
  Mat mat;
  MatStencil row, col[5];
  PetscReal v[3] = {4, -1, -1};
  PetscInt i,j;
  
  PetscFunctionBegin;
  PetscLogEventBegin(EVENT_GenerateLaplacian2DNoBC,0,0,0,0);
  
  ierr = DACreate2d(PETSC_COMM_SELF,//MPI Communicator   
  DA_NONPERIODIC,   //DA_NONPERIODIC, DA_XPERIODIC, DA_YPERIODIC, DA_XYPERIODIC
  DA_STENCIL_STAR, //DA_STENCIL_BOX or DA_STENCIL_STAR
  d1, d2, //Global dimension
  1, 1,          //Number procs per dim
  1,            //dof
  1,           //stencil width
  0,0,        //specific array of nodes
  &da);
  
  DAGetMatrix( da, MATSEQSBAIJ, &mat);
//  MatSetOption(mat, MAT_ROWS_SORTED);
//  MatSetOption(mat, MAT_COLUMNS_SORTED );
    
  for( j = 0; j < d2 - 1; ++j)
  {
    for( i = 0; i < d1 - 1; ++i)
    {
      row.i = i;
      row.j = j;
      col[0].i = i; col[0].j = j;
      col[1].i = i; col[1].j = j+1;
      col[2].i = i+1; col[2].j = j;
      ierr = MatSetValuesStencil(mat, 1, &row, 3, col, v, INSERT_VALUES); CHKERRQ(ierr);
    }
  }
  i = d1 - 1;
  for( j = 0; j < d2 - 1; ++j)
  {
    row.i = i;
    row.j = j;
    col[0].i = i; col[0].j = j;
    col[1].i = i; col[1].j = j+1;
    ierr = MatSetValuesStencil(mat,1,&row,2,col,v,INSERT_VALUES); CHKERRQ(ierr);
  }
  /*
  i = d1 - 1;
  row.i = i;
  col[0].i = i;
  col[1].i = i;  
  for( j = 0; j < d2 - 1; ++j)
  {
    row.j = j;
    col[0].j = j;
    col[1].j = j+1;
    ierr = MatSetValuesStencil(mat,1,&row,2,col,v,INSERT_VALUES); CHKERRQ(ierr);
  }
 */
  j = d2 - 1;
  for( i = 0; i < d1 - 1; ++i)
  {
    row.i = i;
    row.j = j;
    col[0].i = i; col[0].j = j;
    col[1].i = i+1; col[1].j = j;
    ierr = MatSetValuesStencil(mat,1,&row,2,col,v,INSERT_VALUES); CHKERRQ(ierr);
  }
     row.i = d1-1;    row.j = d2-1;
  col[0].i = d1-1; col[0].j = d2-1;
  ierr = MatSetValuesStencil(mat,1,&row,1,col,v,INSERT_VALUES); CHKERRQ(ierr);
  
  MatAssemblyBegin(mat,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(mat,MAT_FINAL_ASSEMBLY);
  *m = mat;
  
  PetscLogEventEnd(EVENT_GenerateLaplacian2DNoBC,0,0,0,0);
  PetscFunctionReturn(0);
}
