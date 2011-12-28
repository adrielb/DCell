#include "Common.h"

int main(int argc, char **args)
{
  PetscErrorCode  ierr;
  ierr = DCellInit(); CHKERRQ(ierr);

  MemCache mc;
  ierr = MemCacheCreate(sizeof(iCoor),10,&mc); CHKERRQ(ierr);

  const int len = 30;
  iCoor *coor[len];

  int i;
  printf("Testing Alloc()\n");
  for ( i = 0; i < len/2; ++i) {
    ierr = MemCacheAlloc(mc,&coor[i]); CHKERRQ(ierr);
    printf("%d [%d,%d]\n",i,coor[i]->x,coor[i]->y);
  }

  printf("Testing Free()\n");
  for ( i = 0; i < len/2; ++i) {
    printf("free(%d)\n",i);
    ierr = MemCacheFree(mc,coor[i]); CHKERRQ(ierr);
  }

  printf("Testing Alloc()\n");
  for ( i = 0; i < len; ++i) {
    ierr = MemCacheAlloc(mc,&coor[i]); CHKERRQ(ierr);
    coor[i]->x = i;
    coor[i]->y = 2*i;
    printf("%d [%d,%d]\n",i,coor[i]->x,coor[i]->y);
  }

  ierr = MemCacheDestroy(mc); CHKERRQ(ierr);

  ierr = DCellFinalize(); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
