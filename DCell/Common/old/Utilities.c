#include "Utilities.h"

struct _WriteVec {
  char *format;
  char *dir;
  int count; // file count
  int len; // vector length
  Vec v;
};

/* TODO: parallel save for a global vector by saving only the local part to local disk
 * DAGlobalToLocal( global, local );
 * VecView( COMM_SELF, local ); 
 * 
 * TODO: then gather to rank zero at then end of the simulation
 * VecLoad( "local_name_on_disk", local); 
 * DALocalToGlobal( local, global);
 * VecView( COMM_WORLD, global);
 */

int WriteVecGetCount(WriteVec wv) {  return wv->count;  }

#undef __FUNCT__
#define __FUNCT__ "WriteVecCreate"
PetscErrorCode WriteVecCreate( Vec v, char* name, WriteVec *write )
{
  WriteVec wv;
  int s = 256;
  char temp_dir[256];
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  ierr = PetscNew(struct _WriteVec, &wv); CHKERRQ(ierr);
  ierr = PetscMalloc(s*sizeof(char), &wv->dir); CHKERRQ(ierr);
  ierr = PetscMalloc(s*sizeof(char), &wv->format); CHKERRQ(ierr);
  
  ierr = PetscGetTmp( PETSC_COMM_WORLD, (char*)&temp_dir, s); CHKERRQ(ierr);
  
  sprintf(wv->format, "%s/%s.%%d.Real64",temp_dir, name);
  
  wv->v = v;
  
  *write = wv;
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "WriteVecDestroy"
PetscErrorCode WriteVecDestroy( WriteVec wv )
{
  PetscErrorCode ierr;
 
  PetscFunctionBegin;
  
  ierr = PetscFree( wv->dir ); CHKERRQ(ierr);
  ierr = PetscFree( wv->format ); CHKERRQ(ierr);
  ierr = PetscFree( wv ); CHKERRQ(ierr);  
  
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "WriteVecToDisk"
PetscInt EVENT_WriteVecToDisk;
PetscErrorCode WriteVecToDisk( WriteVec wv )
{
  PetscReal *data;
  int fd;
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  PetscLogEventBegin(EVENT_WriteVecToDisk,0,0,0,0);
  VecGetArray(wv->v,&data);
  VecGetSize(wv->v, &wv->len);
  sprintf(wv->dir, wv->format, wv->count);
//  printf("\ndir: %s\n", wv->dir);
  ierr = PetscBinaryOpen(wv->dir, FILE_MODE_WRITE, &fd); CHKERRQ(ierr);
  ierr = PetscBinaryWrite(fd,data,wv->len,PETSC_REAL,PETSC_FALSE); CHKERRQ(ierr);
  ierr = PetscBinaryClose(fd); CHKERRQ(ierr);
  wv->count++;
  
  VecRestoreArray(wv->v,&data);
  PetscLogEventEnd(EVENT_WriteVecToDisk,0,0,0,0);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "WriteVector"
PetscInt EVENT_WriteVector;
PetscErrorCode WriteVector( char *name, Vec v )
{
  PetscErrorCode ierr;
  PetscReal *data;
  PetscInt len;
  
  PetscFunctionBegin;
  VecGetArray(v,&data);
  VecGetSize(v,&len);
  
  ierr = WriteVectorArray(name, len, data); CHKERRQ(ierr);
  
  VecRestoreArray(v,&data);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "WriteVectorN"
PetscErrorCode WriteVectorN( char *name, int i, Vec v)
{
  PetscErrorCode ierr;
  PetscReal *data;
  char filename[256];
  PetscInt len;
  
  PetscFunctionBegin;
  VecGetArray(v,&data);
  VecGetSize(v,&len);
  
  sprintf(filename, "%s.%d", name, i);
  ierr = WriteVectorArray(filename, len, data); CHKERRQ(ierr);
  
  VecRestoreArray(v,&data);
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "WriteVectorArray"
PetscInt EVENT_WriteVectorArray;
PetscErrorCode WriteVectorArray(char *name, PetscInt len, PetscReal *data)
{
  PetscErrorCode ierr;
  int fd;
  char tempdir[256], dir[256];
  size_t dir_len = 256*sizeof(char);
  
  PetscFunctionBegin;
  PetscLogEventBegin(EVENT_WriteVectorArray,0,0,0,0);
  
  ierr = PetscGetTmp( PETSC_COMM_WORLD, tempdir, dir_len); CHKERRQ(ierr);
  ierr = PetscSNPrintf(dir, dir_len, "%s/%s.Real64", tempdir, name); CHKERRQ(ierr);
  ierr = PetscBinaryOpen(dir, FILE_MODE_WRITE, &fd); CHKERRQ(ierr);
  ierr = PetscBinaryWrite(fd,data,len,PETSC_REAL,PETSC_FALSE); CHKERRQ(ierr);
  ierr = PetscBinaryClose(fd); CHKERRQ(ierr);
  
  PetscLogEventEnd(EVENT_WriteVectorArray,0,0,0,0);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MatWrite"
PetscErrorCode MatWrite( char *name, Mat mat )
{
  PetscViewer view;
  int dir_len = 512;
  char tempdir[512], dir[512];
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  
  ierr = PetscGetTmp( PETSC_COMM_WORLD, tempdir, dir_len); CHKERRQ(ierr);
  ierr = PetscSNPrintf(dir, dir_len, "%s/%s.mat", tempdir, name); CHKERRQ(ierr);
  ierr = PetscViewerASCIIOpen(PETSC_COMM_SELF, dir, &view); CHKERRQ(ierr);
  ierr = MatView(mat, view); CHKERRQ(ierr);
  ierr = PetscViewerDestroy(view); CHKERRQ(ierr);
  
  PetscFunctionReturn(0);
}


/*
#include "mkl_blas.h"
void LibraryVersions(void)
{
  int len=198;
  char buf[198];
  MKLGetVersionString(buf, len);
  printf("%s\n",buf);
  PetscGetVersion(buf, len);
  printf("%s\n",buf);
  //also get compiler version, mpi and 
  printf("\n");
}*/