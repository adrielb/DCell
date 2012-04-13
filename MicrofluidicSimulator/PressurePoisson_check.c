START_TEST( AddingMultipleRowsToMatrix )
{
  PetscErrorCode  ierr;
  UserContext     *uc;
  ierr = RegisterEvents(); CHKERRQ(ierr);
  ierr = PetscNew( UserContext, &uc ); CHKERRQ(ierr);
  ierr = InterpretOptions( uc ); CHKERRQ(ierr);
  ierr = ReadFile(uc->imageFileName, uc->n, &uc->filedata); CHKERRQ(ierr);
  ierr = DetermineScale( uc ); CHKERRQ(ierr);
  PetscPrintf(PETSC_COMM_SELF,"-DELTA_X\t%f\n", uc->DELTA_X);
  Histogram( uc );
  ierr = IndexFreeNodes( uc ); CHKERRQ(ierr);
  printf("DOF: %d\n", uc->numNodes);
  ierr = InitializeVectors( uc ); CHKERRQ(ierr);
  ierr = AssembleCoupledStokesEq(uc); CHKERRQ(ierr);
}
END_TEST

START_TEST( SolvePressure_test )
{
  PetscErrorCode ierr;
  UserContext uc_stack;
  UserContext *uc = &uc_stack;
  
  ierr = InterpretOptions(uc); CHKERRQ(ierr);
  ierr = ReadFile(uc->imageFileName, uc->n, &uc->filedata); CHKERRQ(ierr);
  ierr = IndexFreeNodes(uc); CHKERRQ(ierr);
  ierr = InitializeVectors( uc ); CHKERRQ(ierr);
  ierr = AssemblePressureMatrx(uc); CHKERRQ(ierr);
  ierr = AssemblePressureRHS( uc ); CHKERRQ(ierr);
  ierr = SolvePressure( uc ); CHKERRQ(ierr);
  PetscReal imageData[uc->n], *sol;
  PetscMemzero(imageData, uc->n*sizeof(PetscReal));
  VecGetArray(uc->p, &sol);
  for (int i = 0; i < uc->numNodes; ++i)
    imageData[uc->nodes[i].imageIndex] = sol[i];
  VecRestoreArray(uc->p, &sol);
  for( int i = 0; i < uc->n; i++)
    if( i%uc->numcols == 0 )
      printf("\n%f",imageData[i]);
    else
      printf("\t%f",imageData[i]);
  printf("\n\n");
//  MatView(uc->A, PETSC_VIEWER_STDOUT_SELF);
//  VecView(uc->b, PETSC_VIEWER_STDOUT_SELF);
}
END_TEST

START_TEST( AssemblePressureRHS_test )
{
  PetscErrorCode ierr;
  UserContext uc_stack;
  UserContext *uc = &uc_stack;
  
  ierr = InterpretOptions(uc); CHKERRQ(ierr);
  ierr = ReadFile(uc->imageFileName, uc->n, &uc->filedata); CHKERRQ(ierr);
  ierr = IndexFreeNodes(uc); CHKERRQ(ierr);
  ierr = InitializeVectors( uc ); CHKERRQ(ierr);
  ierr = AssemblePressureMatrx(uc); CHKERRQ(ierr);
  ierr = AssemblePressureRHS( uc ); CHKERRQ(ierr);
  VecView( uc->b, PETSC_VIEWER_STDOUT_SELF );
}
END_TEST

START_TEST( AssemblePressureMatrix_test )
{
  PetscErrorCode ierr;
  UserContext uc_stack;
  UserContext *uc = &uc_stack;
  
  ierr = InterpretOptions(uc); CHKERRQ(ierr);
  ierr = ReadFile(uc->imageFileName, uc->n, &uc->filedata); CHKERRQ(ierr);
  ierr = IndexFreeNodes(uc); CHKERRQ(ierr);
  ierr = AssemblePressureMatrx(uc); CHKERRQ(ierr);
  
  ierr = PetscViewerSetFormat(PETSC_VIEWER_STDOUT_SELF, PETSC_VIEWER_ASCII_DENSE); CHKERRQ(ierr);
//  ierr = MatView(uc->A, PETSC_VIEWER_STDOUT_SELF); CHKERRQ(ierr);

  Vec ones, temp1, temp2;
  VecCreateSeq(PETSC_COMM_WORLD, uc->numNodes, &ones);
  VecSet(ones, 1);
  VecDuplicate(ones, &temp1);
  VecDuplicate(ones, &temp2);
  
  MatMult( uc->A, ones, temp1);
  MatMultTranspose( uc->A, ones, temp2);
  
  PetscTruth flg;
  VecEqual(temp1, temp2, &flg);
  fail_unless(flg, "Not Symmetric");
  
  PetscScalar t1,t2;
  VecSum(temp1, &t1);
  VecSum(temp2, &t2);
  fail_unless(t1 != 0, "indefinite");
  fail_unless(t2 != 0, "indefinite");
  VecDestroy(ones);
  VecDestroy(temp1);
  VecDestroy(temp2);
}
END_TEST
