#include "FileIO.h"

#undef __FUNCT__
#define __FUNCT__ "ReadReal64"
PetscLogEvent EVENT_ReadReal64;
PetscErrorCode ReadReal64( char *optName, Grid g )
{
  int fd, LEN = 256;
  char name[256];
  iCoor size = g->n;
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  PetscLogEventBegin(EVENT_ReadReal64,0,0,0,0);
//  PetscLogEventRegister("Read", 0, &EVENT_ReadReal64);
  ierr = PetscOptionsGetString(0,optName,name,LEN,0); CHKERRQ(ierr);
  ierr = PetscBinaryOpen(name,FILE_MODE_READ,&fd); CHKERRQ(ierr);
  ierr = PetscBinaryRead(fd, g->v1, size.x*size.y, PETSC_DOUBLE); CHKERRQ(ierr);
  ierr = PetscBinaryClose(fd); CHKERRQ(ierr);
  PetscLogFlops( 0 );
  PetscLogEventEnd(EVENT_ReadReal64,0,0,0,0);
  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "ReadCascade512File"
PetscInt EVENT_ReadCascade512File;
PetscErrorCode ReadCascade512File(char *name, Grid *img)
{
  PetscErrorCode  ierr;
  int             fd;
  short           *data_short;
  int             d = 512;
  int             dd= d*d;
  Grid            g;
  iCoor           s = {d,d,0};
  
  PetscFunctionBegin;
  PetscLogEventBegin(EVENT_ReadCascade512File,0,0,0,0);
  ierr = PetscMalloc( dd * sizeof(short), &data_short); CHKERRQ(ierr);
  ierr = PetscBinaryOpen(name,FILE_MODE_READ,&fd); CHKERRQ(ierr);
  ierr = PetscBinaryRead(fd, data_short, dd, PETSC_SHORT); CHKERRQ(ierr);
  ierr = PetscBinaryClose(fd); CHKERRQ(ierr);
  ierr = GridCreate(s,&g); CHKERRQ(ierr);
  for( int i = 0; i < dd; i++)
    g->v1[i] = (PetscReal)data_short[i];
  ierr = PetscFree(data_short); CHKERRQ(ierr);  
  
  *img = g;
  
  PetscLogEventEnd(EVENT_ReadCascade512File,0,0,0,0);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "WriteReal64"
PetscLogEvent EVENT_WriteReal64;
PetscErrorCode WriteReal64( char *optName, Grid g )
{
  int fd, LEN = 256;
  char name[256];
  iCoor size = g->n;
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  PetscLogEventBegin(EVENT_WriteReal64,0,0,0,0);
//  PetscLogEventRegister("WriteReal64", 0, &EVENT_WriteReal64);
  ierr = PetscOptionsGetString(0,optName,name,LEN,0); CHKERRQ(ierr);
  ierr = PetscBinaryOpen(name,FILE_MODE_WRITE,&fd); CHKERRQ(ierr);
  ierr = PetscBinaryWrite(fd, g->v1, size.x*size.y, PETSC_DOUBLE,PETSC_TRUE); CHKERRQ(ierr);
  ierr = PetscBinaryClose(fd); CHKERRQ(ierr);
  PetscLogFlops( 0 );
  PetscLogEventEnd(EVENT_WriteReal64,0,0,0,0);
  PetscFunctionReturn(0);
}