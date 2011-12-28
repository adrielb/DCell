#include "Main.h"
#include "ActiveContoursWithoutEdges.h"
   
#undef __FUNCT__
#define __FUNCT__ "PetscMain"
PetscInt EVENT_PetscMain;
PetscErrorCode PetscMain()
{
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  PetscLogEventRegister(&EVENT_PetscMain,"PetscMain", 0);
  PetscLogEventBegin(EVENT_PetscMain,0,0,0,0);
  PetscErrorCode RegisterEvents_ACWE();
  
  
  char *file_format = "/home/share/Images/Joanne/DIR.IGF1 Adaptation-10.spl/image.a.%d.gray";
//  char *file_format = "/home/abergman/Images/image.a.%d.gray";
  char *file_name;
  int file_name_len = 256;
  PetscMalloc( file_name_len * sizeof(char), &file_name);
  PetscSNPrintf(file_name,file_name_len,file_format, 1);
  printf("%s\n", file_name );
  
  Grid2D img;
  
  ierr = ReadCascade512File(file_name, &img); CHKERRQ(ierr);
   
  UserContext *uc;
  
  ierr = UserContextCreate( img, &uc); CHKERRQ(ierr);
  
  WriteVec wv;
  ierr = WriteVecCreate(uc->ls->g2d->v, "phi", &wv); CHKERRQ(ierr);
  
//  ierr = LevelSetInitializeToCircle(uc->ls, 350, 217, 64); CHKERRQ(ierr);
  ierr = Border(img); CHKERRQ(ierr);
  ierr = ThresholdImage(img, 10000.5, uc->ls->g2d); CHKERRQ(ierr);
  ierr = UpdateIrregularNodeList(uc->ls); CHKERRQ(ierr);
  ierr = IrregularNodeListWrite(uc->ls,0); CHKERRQ(ierr);
  ierr = ReinitializeLevelSet(uc->ls); CHKERRQ(ierr);
  
  for( int i = 1; i < 240; i++ )
  {
    printf("i: %d\n", i);
    if( i%10==0 )
    {
      ierr = UpdateIrregularNodeList(uc->ls); CHKERRQ(ierr);
      ierr = IrregularNodeListWrite(uc->ls, WriteVecGetCount(wv) ); CHKERRQ(ierr);
      ierr = ReinitializeLevelSet(uc->ls); CHKERRQ(ierr);
      ierr = WriteVecToDisk(wv); CHKERRQ(ierr);
    }
    
    ierr = SingleStep( uc, img); CHKERRQ(ierr);
  }
  
  PetscLogEventEnd(EVENT_PetscMain,0,0,0,0);
  PetscFunctionReturn(0);
}