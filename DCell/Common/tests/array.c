#include "Common.h"

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc, char **args)
{
  PetscErrorCode  ierr;
  PetscFunctionBegin;
  ierr = PetscInitialize(&argc, &args, (char *) 0, ""); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "Start %s\n", __FILE__); CHKERRQ(ierr);
  
  Array a;
  ierr = ArrayCreate( "test", sizeof(iCoor), 10, &a); CHKERRQ(ierr);
  ierr = ArraySetSize(a, 5); CHKERRQ(ierr);
  
  int i;
  iCoor *c;
  printf("Testing ArrayGet()\n");
  for( i = 0; i < ArrayLength(a); i++ )
  {
    ierr = ArrayGet(a,i,&c); CHKERRQ(ierr);
    c->x = i;
    c->y = 2*i;
  }
  for( i = 0; i < ArrayLength(a); i++ )
  {
    ierr = ArrayGet(a,i,&c); CHKERRQ(ierr);
    printf("%d, %d\n", c->x, c->y);
  }

  printf("Testing ArraySetSize()\n");
  ierr = ArraySetSize(a, 15); CHKERRQ(ierr);
  for( i = 0; i < ArrayLength(a); i++ )
  {
    ierr = ArrayGet(a,i,&c); CHKERRQ(ierr);
    c->x = i;
    c->y = 2*i;
  }
  for( i = 0; i < ArrayLength(a); i++ )
  {
    ierr = ArrayGet(a,i,&c); CHKERRQ(ierr);
    printf("%d, %d\n", c->x, c->y);
  }

  printf("Testing ArrayAppend()\n");
  for( i = 0; i < 16; i++ )
  {
    ierr = ArrayAppend(a,&c); CHKERRQ(ierr);
    c->x = i;
    c->y = 2*i;
  }
  for( i = 0; i < ArrayLength(a); i++ )
  {
    ierr = ArrayGet(a,i,&c); CHKERRQ(ierr);
    printf("%d, %d\n", c->x, c->y);
  }

  printf("Testing Delete1()\n");
  for( i = 0; i < ArrayLength(a); i+=5 )
  {
    ierr = ArrayDelete1(a,i); CHKERRQ(ierr);
  }
  for( i = 0; i < ArrayLength(a); i++ )
  {
    ierr = ArrayGet(a,i,&c); CHKERRQ(ierr);
    printf("%d: %d, %d\n", i, c->x, c->y);
  }

  printf("\nZeroing Array\n");
  ierr = ArrayZero(a); CHKERRQ(ierr);
  for( i = 0; i < ArrayLength(a); i++ )
  {
    ierr = ArrayGet(a,i,&c); CHKERRQ(ierr);
    printf("%d, %d\n", c->x, c->y);
  }

  printf("Destroying Array\n");
  ierr = ArrayDestroy( a ); CHKERRQ(ierr);

  double d = 3.14;
  typedef double* LS;
  LS ls, *lsp;
//  *ls = &d;
  printf("1: %p\t%p\t%p\n", ls, &ls, &d);
  ierr = ArrayCreate( "test2", sizeof(LS), 10, &a); CHKERRQ(ierr);
  ierr = ArraySetSize(a,2); CHKERRQ(ierr);
  ierr = ArrayGet(a,0,(void*)&lsp); CHKERRQ(ierr);
  printf("2: %p\t%p\t%p\n", lsp, &lsp, &d);
  *lsp = &d;
  printf("3: %p\t%p\t%p\n", ls, &ls, &d);
  ierr = ArrayGetP(a,0,(void*)&ls); CHKERRQ(ierr);
  printf("4: %p\t%p\t%p\n", ls, &ls, &d);

  PetscLogDouble usage1,start;
  ierr = PetscMemoryGetCurrentUsage(&usage1); CHKERRQ(ierr);
  PetscReal mb = 100;
  ierr = PetscOptionsGetReal(0,"-mb",&mb,0); CHKERRQ(ierr);
  printf("\nAllocation test: %5.1f MB \n", mb);
  int len = mb * 1048576 / sizeof(double);
  ierr = PetscGetTime(&start); CHKERRQ(ierr);
  ierr = ArrayCreate( "test", sizeof(double), len, &a); CHKERRQ(ierr);
  PetscLogDouble usage2,end;
  ierr = PetscGetTime(&end); CHKERRQ(ierr);
  ierr = PetscMemoryGetCurrentUsage(&usage2); CHKERRQ(ierr);
  printf("U1: %5.1f MB\n", usage1/1048576.);
  printf("U2: %5.1f MB\n", usage2/1048576.);
  printf("time: %f sec\n", end-start);
  ierr = ArrayDestroy( a ); CHKERRQ(ierr);

  printf("\nTest Append Pointer\n");
  ierr = ArrayCreate("testptr",sizeof(void*),10,&a); CHKERRQ(ierr);
  void *ptr1 = (void*)1, *ptr2;
  ierr = ArrayAppendPtr(a,ptr1); CHKERRQ(ierr);
  ierr = ArrayGetP(a,0,&ptr2); CHKERRQ(ierr);
  printf("ptr1: %p\n", ptr1);
  printf("ptr2: %p\n", ptr2);
  ierr = ArrayDestroy(a); CHKERRQ(ierr);

  ierr = PetscPrintf(PETSC_COMM_WORLD, "End %s\n", __FILE__); CHKERRQ(ierr);
  ierr = PetscFinalize(); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
