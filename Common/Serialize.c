#include "Common.h"
#include "Grid.h"

static const char WORKING_DIRECTORY[PETSC_MAX_PATH_LEN];
PetscLogEvent EVENT_GridWrite;

#undef __FUNCT__
#define __FUNCT__ "DCellGetWorkingDirectory"
PetscErrorCode DCellGetWorkingDirectory( char* workingdirectory )
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  if( WORKING_DIRECTORY[0] == 0 ) {
    ierr = PetscGetTmp( PETSC_COMM_WORLD, (char*)WORKING_DIRECTORY, PETSC_MAX_PATH_LEN); CHKERRQ(ierr);
  }
  ierr = PetscStrcpy(workingdirectory,WORKING_DIRECTORY); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "GridWrite"
PetscErrorCode GridWrite( Grid g, int t )
{
  size_t len = PETSC_MAX_PATH_LEN;
  char filename[PETSC_MAX_PATH_LEN];
  char wd[PETSC_MAX_PATH_LEN];
  int time = t + FILE_COUNT_START;
  PetscViewer binv;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscLogEventBegin(EVENT_GridWrite,0,0,0,0); CHKERRQ(ierr);
  ierr = DCellGetWorkingDirectory(wd); CHKERRQ(ierr);
  ierr = PetscSNPrintf(filename,len,"%s/%s.%d.Real64",wd,g->name,time); CHKERRQ(ierr);
  ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF,filename,FILE_MODE_WRITE,&binv); CHKERRQ(ierr);
  ierr = VecView(g->v, binv); CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&binv); CHKERRQ(ierr);

  if( g->filePos == NULL ) {
    ierr = PetscSNPrintf(filename,len,"%s/%s.pos",wd,g->name); CHKERRQ(ierr);
    ierr = PetscViewerASCIIOpen(PETSC_COMM_SELF,filename,&g->filePos); CHKERRQ(ierr);
    ierr = PetscSNPrintf(filename,len,"%s/%s.size",wd,g->name); CHKERRQ(ierr);
    ierr = PetscViewerASCIIOpen(PETSC_COMM_SELF,filename,&g->fileSize); CHKERRQ(ierr);
  }
  if( g->is2D ) {
    ierr = PetscViewerASCIIPrintf(g->filePos,"%d %d\n",g->p.x,g->p.y); CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(g->fileSize,"%d %d\n",g->n.x,g->n.y); CHKERRQ(ierr);
  } else {
    ierr = PetscViewerASCIIPrintf(g->filePos, "%d %d %d\n",g->p.x,g->p.y,g->p.z); CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(g->fileSize,"%d %d %d\n",g->n.x,g->n.y,g->n.z); CHKERRQ(ierr);
  }
  ierr = PetscLogEventEnd(EVENT_GridWrite,0,0,0,0); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "VecWrite"
PetscErrorCode VecWrite( Vec vec, const char *name, int t )
{
  size_t len = PETSC_MAX_PATH_LEN;
  char filename[PETSC_MAX_PATH_LEN];
  char wd[PETSC_MAX_PATH_LEN];
  int time = t + FILE_COUNT_START;
  PetscViewer binv;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = DCellGetWorkingDirectory(wd); CHKERRQ(ierr);
  ierr = PetscSNPrintf(filename,len,"%s/%s.%d.Real64",wd,name,time); CHKERRQ(ierr);
  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,filename,FILE_MODE_WRITE,&binv); CHKERRQ(ierr);
  ierr = VecView(vec, binv); CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&binv); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MatWrite"
PetscErrorCode MatWrite( Mat mat, const char *name, int t )
{
  size_t len = PETSC_MAX_PATH_LEN;
  char tempdir[PETSC_MAX_PATH_LEN];
  char filename[PETSC_MAX_PATH_LEN];
  int time = t + FILE_COUNT_START;
  PetscViewer view;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscGetTmp( PETSC_COMM_WORLD, tempdir, len); CHKERRQ(ierr);
  ierr = PetscSNPrintf(filename,len,"%s/%s.%d.mat",tempdir,name,time); CHKERRQ(ierr);
  ierr = PetscViewerASCIIOpen(PETSC_COMM_SELF, filename, &view); CHKERRQ(ierr);
  ierr = MatView(mat, view); CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&view); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

