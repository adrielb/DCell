#include "Common.h"

PetscErrorCode  RegisterEvents(void);

// source file = __FILE__
const static char MAINFILE[PETSC_MAX_PATH_LEN];

#undef __FUNCT__
#define __FUNCT__ "DCellInitialize"
PetscErrorCode  DCellInitialize(int *argc,char ***args, const char file[])
{
  setbuf(stdout, NULL);
  int stack = 10e6; //TODO: make GA stack/heap a petsc option -ga_stack -ga_heap
  int heap =  100e6;
  PetscErrorCode  ierr;
  ierr = MPI_Init(argc,args); CHKERRQ(ierr);
  if( !MA_init(C_DBL, stack, heap)) {
    GA_Error((char*)"MA_init failed",stack+heap);
  }
  GA_Initialize();
  ierr = PetscInitialize(argc, args, (char *) 0, ""); CHKERRQ(ierr);

  PetscFunctionBegin;
  PetscViewer viewer;
  ierr = PetscViewerCreate(PETSC_COMM_SELF,&viewer); CHKERRQ(ierr);
  ierr = PetscViewerSetType(viewer, PETSCVIEWERASCII); CHKERRQ(ierr);
  ierr = PetscViewerFileSetMode(viewer, FILE_MODE_WRITE); CHKERRQ(ierr);
  ierr = PetscViewerFileSetName(viewer, "params.txt"); CHKERRQ(ierr);
  ierr = PetscOptionsView(viewer); CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);
  //TODO: Open HDF5 simulation file here
  ierr = RegisterEvents(); CHKERRQ(ierr);
  ierr = HeapRegisterEvents(); CHKERRQ(ierr);
  ierr = PetscStrcpy((char*)MAINFILE,file); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "Start %s\n", MAINFILE); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DCellFinalize"
PetscErrorCode  DCellFinalize()
{
  PetscErrorCode  ierr;
  PetscFunctionBegin;

  //TODO: Close HDF5 simulation file here
  GA_Terminate();
  ierr = PetscPrintf(PETSC_COMM_WORLD, "End %s\n", MAINFILE); CHKERRQ(ierr);
  ierr = PetscFinalize(); CHKERRQ(ierr);
  ierr = MPI_Finalize(); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscLogEvent EVENT_GridWrite;
PetscLogEvent EVENT_ArraySetSize;
PetscLogEvent EVENT_GAPutVec;
PetscLogEvent EVENT_GAGetVec;
PetscLogEvent EVENT_GAGather;
PetscLogEvent EVENT_GAScatterAcc;

#undef __FUNCT__
#define __FUNCT__ "RegisterEvents"
PetscErrorCode  RegisterEvents(void)
{
  PetscErrorCode  ierr;
  PetscFunctionBegin;

  ierr = PetscLogEventRegister("GridWrite", 0, &EVENT_GridWrite); CHKERRQ(ierr);
  ierr = PetscLogEventRegister("ArrayResize", 0, &EVENT_ArraySetSize); CHKERRQ(ierr);

  ierr = PetscLogEventRegister("GAPutVec", 0, &EVENT_GAPutVec); CHKERRQ(ierr);
  ierr = PetscLogEventRegister("GAGetVec", 0, &EVENT_GAGetVec); CHKERRQ(ierr);
  ierr = PetscLogEventRegister("GAGather", 0, &EVENT_GAGather); CHKERRQ(ierr);
  ierr = PetscLogEventRegister("GAScatterAcc", 0, &EVENT_GAScatterAcc); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

