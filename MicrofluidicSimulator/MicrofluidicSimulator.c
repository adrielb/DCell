#include "MicrofluidicSimulator.h"

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc, char **args)
{
  PetscErrorCode  ierr;
  PetscFunctionBegin;
  ierr = PetscInitialize(&argc, &args, (char *) 0, ""); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "Start %s\n", __FILE__); CHKERRQ(ierr);

  ierr = PetscPrintf(PETSC_COMM_WORLD, "End %s\n", __FILE__); CHKERRQ(ierr);
  ierr = PetscFinalize(); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PetscMain"
PetscErrorCode PetscMain()
{
  PetscErrorCode  ierr;
  UserContext     *uc;
  
  MicrofluidicSimulator_RegisterEvents();
  
  ierr = PetscNew( UserContext, &uc ); CHKERRQ(ierr);
  ierr = InterpretOptions( uc ); CHKERRQ(ierr);
  ierr = ReadFile(uc->imageFileName, uc->n, &uc->filedata); CHKERRQ(ierr);
  ierr = DetermineScale( uc ); CHKERRQ(ierr);
  Histogram( uc );
  ierr = IndexFreeNodes( uc ); CHKERRQ(ierr);
  
  /*  Determine file info  */
  printf("-n\t\t\t%d\n", uc->n);
  printf("DOF: %d\n", uc->numNodes);
  
  Vec u,v;
  VecCreateSeq(PETSC_COMM_WORLD, uc->numNodes, &u);
  VecDuplicate(u, &v);
  
//  PressureIncrement(uc, u, v);
  
  WriteResult(uc,u,"u_vel.Real64");
  WriteResult(uc,v,"v_vel.Real64");
  
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DetermineScale"
PetscInt EVENT_DetermineScale;
PetscErrorCode DetermineScale( UserContext *uc )
{
  PetscErrorCode ierr;
  int i, j;
  PetscInt maxX=0, maxY=0, minX=uc->numcols, minY=uc->numrows;
  
  PetscFunctionBegin;
  PetscLogEventBegin(EVENT_DetermineScale,0,0,0,0);
  for( j = 0; j < uc->numrows; ++j )
  {
    for( i = 0; i < uc->numcols; ++i)
    { 
      if( uc->filedata[i + j*uc->numcols] == uc->SCALE_COLOR )
      {
        maxX = MAX(i, maxX);
        maxY = MAX(j, maxY);
        minX = MIN(i, minX);
        minY = MIN(j, minY);
      }
    }
  }
  
  uc->DELTA_X = uc->SCALE_BAR / ( (double)MAX(maxX - minX, maxY - minY) );
  
  if( uc->DELTA_X < 0. ) uc->DELTA_X = 1;

  printf("-SCALE_COLOR:\t%d\n", uc->SCALE_COLOR);
  printf("-SCALE_BAR:\t%f\n", uc->SCALE_BAR);
  printf("maxX: %d\n",maxX);
  printf("minX: %d\n",minX);
  printf("maxY: %d\n",maxY);
  printf("minY: %d\n",minY);
  printf("dX:\t%d\n", maxX-minX);
  printf("dY:\t%d\n", maxY-minY);
  printf("-DELTA_X\t%f\n", uc->DELTA_X);
  
  PetscLogEventEnd(EVENT_DetermineScale,0,0,0,0);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "Histogram"
PetscInt EVENT_Histogram;
PetscErrorCode Histogram( UserContext *uc )
{
  int hist[256];
  PetscErrorCode ierr;
  
  
  PetscFunctionBegin;
  PetscLogEventBegin(EVENT_Histogram,0,0,0,0);
  PetscLogEventRegister(&EVENT_Histogram,"Histogram", 0);
  
  for( int i = 0; i < 256; i++)
    hist[i] = 0;
  
  for( int i = 0; i < uc->n; i++)
    hist[uc->filedata[i]]++;
  
  printf("\tCOLOR\tCOUNT\n");
  for( int i = 0; i < 256; i++)
  {
    if( hist[i] == 0 ) continue;
    printf("\t%d\t\t%d\n", i, hist[i]);
  }
  printf("\n");
  
  PetscLogEventEnd(EVENT_Histogram,0,0,0,0);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "InterpretOptions"
PetscErrorCode InterpretOptions(UserContext *uc)
{
  PetscErrorCode ierr;
  PetscTruth flg;
  
  PetscFunctionBegin;
  PetscMalloc(256*sizeof(char), &uc->imageFileName);
  PetscOptionsHasName( PETSC_NULL, "-check", &flg);
  if( flg )
  {
    ierr = PetscOptionsInsertFile(PETSC_COMM_WORLD,"check.params",PETSC_FALSE); CHKERRQ(ierr);
  } else {
    ierr = PetscOptionsInsertFile(PETSC_COMM_WORLD,"StokesFlow.params",PETSC_FALSE); CHKERRQ(ierr);
  }
  ierr = PetscOptionsGetString(PETSC_NULL, "-IMAGE_NAME",           uc->imageFileName, 256, 0); CHKERRQ(ierr);
  ierr = PetscOptionsGetInt(   PETSC_NULL, "-HEIGHT",              &uc->numrows, 0); CHKERRQ(ierr);
  ierr = PetscOptionsGetInt(   PETSC_NULL, "-WIDTH",               &uc->numcols, 0); CHKERRQ(ierr);
  ierr = PetscOptionsGetInt(   PETSC_NULL, "-BACKGROUND_COLOR",    &uc->BACKGROUND_COLOR, 0); CHKERRQ(ierr);
  ierr = PetscOptionsGetInt(   PETSC_NULL, "-FLUIDIC_LAYER_COLOR", &uc->FLUIDIC_LAYER_COLOR, 0); CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(  PETSC_NULL, "-VISCOSITY",           &uc->VISCOSITY, 0); CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(  PETSC_NULL, "-DIFFUSION",           &uc->DIFFUSION, 0); CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(  PETSC_NULL, "-DELTA_T",             &uc->DELTA_T, 0); CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(  PETSC_NULL, "-CHIP_HEIGHT",         &uc->DELTA_Z, 0); CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(  PETSC_NULL, "-SCALE_BAR",           &uc->SCALE_BAR, 0); CHKERRQ(ierr);
  ierr = PetscOptionsGetInt(   PETSC_NULL, "-SCALE_COLOR",         &uc->SCALE_COLOR, 0); CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(  PETSC_NULL, "-DENSITY",             &uc->DENSITY, 0); CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(  PETSC_NULL, "-GRAVITY",             &uc->GRAVITY, 0); CHKERRQ(ierr);  

  char buf[10];
  for( int i = 0; i < 256; i++)
  {
    sprintf(buf, "-%d", i);
    ierr = PetscOptionsGetReal("BC_HEIGHT_",        buf, &uc->BCLabels[i].pressureBC,0); CHKERRQ(ierr);
    ierr = PetscOptionsGetReal("BC_CONCENTRATION_", buf, &uc->BCLabels[i].concentrationBC,0); CHKERRQ(ierr);
    uc->BCLabels[i].pressureBC *= uc->DENSITY * uc->GRAVITY;
  }
  
  uc->n = uc->numrows * uc->numcols;
  uc->DELTA_Z /= 2;

  ierr = PetscMalloc(uc->n*sizeof(PetscReal),&uc->imageResult); CHKERRQ(ierr);
  
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "IndexFreeNodes"
PetscInt EVENT_IndexFreeNodes;
PetscErrorCode IndexFreeNodes(UserContext* uc)
{
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  PetscLogEventBegin(EVENT_IndexFreeNodes, 0,0,0,0);
//Index all nodes, {0, ..., numNodes}
  uc->numNodes = 0;
  ierr = PetscMalloc(uc->n*sizeof(PetscInt),&uc->imageToNode); CHKERRQ(ierr);
  for (int i = 0; i < uc->n; ++i)
  {
    if( uc->filedata[i] != uc->BACKGROUND_COLOR && 
        uc->filedata[i] != uc->SCALE_COLOR ) // Either FLUIDIC or BC
    {
      uc->imageToNode[i] = uc->numNodes;
      uc->numNodes++;
    } else {
      uc->imageToNode[i] = -1;
    }
  }
  
  Node *n;
  ierr = PetscMalloc(uc->numNodes*sizeof(Node),&uc->nodes); CHKERRQ(ierr);
  PetscInt col[4], i , j, count = 0, nodeIndex;
  const PetscInt C[4] = {-uc->numcols, -1, 1, uc->numcols};
  const PetscInt BOX[9] = {
    -uc->numcols-1, -uc->numcols, -uc->numcols+1,
                -1,            0,              1,
     uc->numcols-1,  uc->numcols,  uc->numcols+1};
  
  for( i = 0; i < uc->n; ++i) 
  {
    if( uc->imageToNode[i] == -1 )  continue;
    n = &uc->nodes[count];
    n->imageIndex = i;
    n->numNei = 0;
    n->isInterior = PETSC_TRUE;
        
    for( j = 0; j < 4; ++j)
    {
      col[j] = i + C[j];
      n->star[j] = -1;
      if(uc->filedata[col[j]] != uc->BACKGROUND_COLOR )
      {
        nodeIndex = uc->imageToNode[col[j]];
        n->star[j] = nodeIndex;
        n->nei[n->numNei] = nodeIndex;
        n->numNei++;
      }
    }
    for( j = 0; j < 9; j++)
    {
      if( uc->filedata[ i+BOX[j] ] == uc->BACKGROUND_COLOR )
      {
        n->isInterior = PETSC_FALSE;
        break;
      }
    }
    n->nei[n->numNei] = count;
    n->star[4] = count;
    count++;
  }
  
//Count the number of BC nodes
  uc->numBC = 0;
  for( int i = 0; i < uc->n; i++ )
  {
    if( uc->filedata[i] != uc->BACKGROUND_COLOR && 
        uc->filedata[i] != uc->FLUIDIC_LAYER_COLOR &&
        uc->filedata[i] != uc->SCALE_COLOR )    
        uc->numBC++;
  }
  ierr = PetscMalloc(uc->numBC*sizeof(BCNode), &uc->bcNodes); CHKERRQ(ierr);

//Index the BC nodes to their nodal index
  count = 0;
  for( int i = 0; i < uc->n; i++ )
  {
    if( uc->filedata[i] != uc->BACKGROUND_COLOR && 
        uc->filedata[i] != uc->FLUIDIC_LAYER_COLOR &&
        uc->filedata[i] != uc->SCALE_COLOR )
    {    
      uc->bcNodes[count].nodeIndex = uc->imageToNode[i];
      uc->bcNodes[count].pressureBC = uc->BCLabels[ uc->filedata[i] ].pressureBC;
      uc->bcNodes[count].concentrationBC = uc->BCLabels[ uc->filedata[i] ].concentrationBC;
      count++;
    } 
  }
  PetscLogEventEnd(EVENT_IndexFreeNodes, 0,0,0,0);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "WriteResult"
PetscInt EVENT_WriteResult;
PetscErrorCode WriteResult( UserContext *uc, Vec v, char *filename )
{
  PetscErrorCode ierr;
  PetscReal      *sol;
  int            fd;
  
  PetscFunctionBegin;
  PetscLogEventBegin(EVENT_WriteResult,0,0,0,0);
  ierr = PetscMemzero(uc->imageResult, uc->n * sizeof(PetscReal)); CHKERRQ(ierr);
  VecGetArray(v, &sol);
  
  for (int i = 0; i < uc->numNodes; ++i)
    uc->imageResult[uc->nodes[i].imageIndex] = sol[i];
  
//  ierr = WriteVectorArray( filename, uc->n, uc->imageResult); CHKERRQ(ierr);

  VecRestoreArray(v,&sol);
  PetscLogEventEnd(EVENT_WriteResult,0,0,0,0);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DestroyContext"
PetscErrorCode DestroyContext( UserContext* uc )
{
  PetscErrorCode ierr;
  
  PetscFunctionBegin;

  ierr = PetscFree(uc->imageFileName); CHKERRQ(ierr);
  ierr = PetscFree(uc->filedata); CHKERRQ(ierr);
  ierr = PetscFree(uc->imageToNode); CHKERRQ(ierr);
  ierr = PetscFree(uc->nodes); CHKERRQ(ierr);
  ierr = PetscFree(uc->imageResult); CHKERRQ(ierr);
  ierr = PetscFree(uc); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ReadFile"
PetscInt EVENT_ReadFile;
PetscErrorCode ReadFile(char *name, PetscInt len, unsigned char **filedata)
{
  PetscErrorCode ierr;
  int fd;
  unsigned char *data;
  
  PetscFunctionBegin;
  PetscLogEventBegin(EVENT_ReadFile,0,0,0,0);
  
  ierr = PetscMalloc( len, &data); CHKERRQ(ierr);
  ierr = PetscBinaryOpen(name,FILE_MODE_READ,&fd); CHKERRQ(ierr);
  ierr = PetscBinaryRead(fd, data, len, PETSC_CHAR); CHKERRQ(ierr);
  ierr = PetscBinaryClose(fd); CHKERRQ(ierr);
  
  *filedata = data;
 
  PetscLogEventEnd(EVENT_ReadFile,0,0,0,0);
  PetscFunctionReturn(0);
}

void MicrofluidicSimulator_RegisterEvents()
{
//  RegisterUtilityEvents();
  PetscLogEventRegister(&EVENT_ReadFile," ReadFile", 0);
  PetscLogEventRegister(&EVENT_WriteResult," WriteResult", 0);
  PetscLogEventRegister(&EVENT_IndexFreeNodes, " IndexFreeNodes", 0);
//  PetscLogEventRegister(&EVENT_SolvePressure," SolvePressure", 0);  
//  PetscLogEventRegister(&EVENT_ComputeShearStress," ShearStress", 0);
  PetscLogEventRegister(&EVENT_DetermineScale," DetermineScale", 0);
//  PetscLogEventRegister(&EVENT_AssembleDiffusionMatrix,"AssembleDiffusionMatrix", 0);  
//  PetscLogEventRegister(&EVENT_AssembleConcentrationRHS,"AssembleConcentrationRHS", 0);
//  PetscLogEventRegister(&EVENT_SolveConcentration,"SolveConcentration", 0);
}
