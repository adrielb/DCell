#include "Common.h"

struct _Array {
  MPI_Comm comm;
  void *dataArray;
  int len;        // actual length of data in array
  int ELEMSIZE;   // a constant
  size_t MAXSIZE; // current max allocate number of elements
  PetscReal scale;// scaling factor when resizing array beyond requested size
  char name[64];  // name used for PetscInfo()
  iCoor p,q;
  iCoor size;
};

PetscLogEvent EVENT_ArraySetSize;

#undef __FUNCT__
#define __FUNCT__ "ArrayCreate"
PetscErrorCode ArrayCreate( const char name[], int elemSize, Array *array )
{
  Array a;
  int initSize = 1024;
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  ierr = PetscNew(struct _Array, &a); CHKERRQ(ierr);
  ierr = PetscOptionsGetInt( name, "-array_initsize", &initSize, 0); CHKERRQ(ierr);
  ierr = PetscMalloc( elemSize*initSize, &a->dataArray); CHKERRQ(ierr);
  ierr = PetscMemzero(a->dataArray, elemSize*initSize); CHKERRQ(ierr);
  a->ELEMSIZE = elemSize;
  a->MAXSIZE = initSize;
  a->len = 0;
  a->scale = 1.1;
  ierr = PetscStrcpy(a->name, name); CHKERRQ(ierr);
  ierr = PetscOptionsGetReal( name, "-array_scale", &a->scale, 0); CHKERRQ(ierr);
  ierr = PetscInfo2(0, "%s initially %d MB\n", name, elemSize*initSize / (1024*1024)); CHKERRQ(ierr);
  
  *array = a;
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ArrayDestroy"
PetscErrorCode ArrayDestroy( Array a )
{
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  if( a->dataArray )
  {
    ierr = PetscFree( a->dataArray ); CHKERRQ(ierr);
  }

  ierr = PetscFree( a ); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ArrayDuplicate"
PetscErrorCode ArrayDuplicate( Array a, Array *newarray )
{
  Array new;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = ArrayCreate(a->name, a->ELEMSIZE, &new); CHKERRQ(ierr);
  ierr = ArraySetSize(new, a->MAXSIZE ); CHKERRQ(ierr);
  *newarray = new;
  PetscFunctionReturn(0);
}

//TODO: Time resizing events, frequency, and size
#undef __FUNCT__
#define __FUNCT__ "ArraySetSize"
PetscErrorCode ArraySetSize( Array a, int size )
{
  PetscErrorCode ierr;
  
  PetscFunctionBegin;

  if( size < 0 )
  {
    SETERRQ2(a->comm, PETSC_ERR_ARG_OUTOFRANGE,"Array [%s] size %d < 0", a->name, size);
  }

  a->len = size;
  
  if( a->MAXSIZE < size ) //TODO: report resizing in petsc info
  {
    ierr = PetscLogEventBegin(EVENT_ArraySetSize,0,0,0,0); CHKERRQ(ierr);
    int s = a->scale * size * a->ELEMSIZE;
    void *tmp;
    ierr = PetscInfo4(0,"%s resizing: %d to %d (%d MB)\n",a->name, a->MAXSIZE, (int)(a->scale*size), s/(1024*1024) ); CHKERRQ(ierr);
    ierr = PetscMalloc( s, &tmp ); CHKERRQ(ierr);
    ierr = PetscMemzero(tmp, s); CHKERRQ(ierr); //TODO: is this redundant?
    ierr = PetscMemcpy(tmp,a->dataArray,a->ELEMSIZE*a->MAXSIZE); CHKERRQ(ierr);
    ierr = PetscFree(a->dataArray); CHKERRQ(ierr);
    a->dataArray = tmp;
    a->MAXSIZE = s / a->ELEMSIZE;
    ierr = PetscLogEventEnd(EVENT_ArraySetSize,0,0,0,0); CHKERRQ(ierr);
  }
  
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ArrayWrite"
PetscErrorCode ArrayWrite( Array a, int t )
{
  size_t len = 512;
  char tempdir[512];
  char filename[512];
  int fp;
  int time = t + FILE_COUNT_START;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscGetTmp( PETSC_COMM_WORLD, tempdir, len); CHKERRQ(ierr);
  ierr = PetscSNPrintf(filename,len,"%s/%s.%d.array",tempdir,a->name,time); CHKERRQ(ierr);
  ierr = PetscBinaryOpen(filename,FILE_MODE_WRITE,&fp); CHKERRQ(ierr);
  ierr = PetscBinaryWrite(fp,a->dataArray,a->len*a->ELEMSIZE,PETSC_CHAR,PETSC_FALSE); CHKERRQ(ierr);
  ierr = PetscBinaryClose(fp); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ArrayAppend"
PetscErrorCode ArrayAppend( Array a, void *elem )
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = ArraySetSize( a, a->len+1); CHKERRQ(ierr);
  ierr = ArrayGet( a, a->len-1,(void**)elem ); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ArrayAppendPtr"
PetscErrorCode ArrayAppendPtr( Array a, void *elem ) // append input pointer to array of pointers
{
  void **last;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = ArraySetSize( a, a->len+1); CHKERRQ(ierr);
  ierr = ArrayGet( a, a->len-1,&last ); CHKERRQ(ierr);
  *last = elem;
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ArrayGet"
PetscErrorCode ArrayGet( Array a, int i, void *elem )
{
  PetscFunctionBegin;
  if( 0 <= i && i < a->len )
    *((void**)elem) = a->dataArray + a->ELEMSIZE * i;
  else
    SETERRQ3(a->comm, PETSC_ERR_ARG_OUTOFRANGE,"ArrayGet[%s]: Index %d not in [0 - %d)",a->name,i,a->len);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ArrayGetP"
PetscErrorCode ArrayGetP( Array a, int i, void *elem )
{
  PetscFunctionBegin;
  if( 0 <= i && i < a->len )
    *(void**)elem = *(void**)(a->dataArray + a->ELEMSIZE * i);
  else
    SETERRQ3(a->comm, PETSC_ERR_ARG_OUTOFRANGE,"ArrayGet[%s]: Index %d not in [0 - %d)",a->name,i,a->len);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ArrayGetCoor"
PetscErrorCode ArrayGetCoor( Array a, iCoor pos, void *elem)
{
  int idx;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  idx = (pos.x - a->p.x) + a->size.x * (pos.y - a->p.y) + a->size.x * a->size.y * (pos.z - a->p.z);
  ierr = ArrayGet(a, idx, elem); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ArrayGetCoorP"
PetscErrorCode ArrayGetCoorP( Array a, iCoor pos, void *elem)
{
  int idx;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  idx = (pos.x - a->p.x) + a->size.x * (pos.y - a->p.y) + a->size.x * a->size.y * (pos.z - a->p.z);
  ierr = ArrayGetP(a, idx, elem); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ArraySetCoor"
PetscErrorCode ArraySetCoor( Array a, iCoor shift, iCoor size )
{
  PetscErrorCode ierr;
  PetscFunctionBegin;
  if( size.z == 0 ) size.z = 1;
  a->p = shift;
  a->size = size;
  a->q.x = shift.x + size.x;
  a->q.y = shift.y + size.y;
  a->q.z = shift.z + size.z;
  ierr = ArraySetSize(a,size.x*size.y*size.z); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

int ArrayLength( Array a )
{
  return a->len;
}

size_t ArrayMaxSize( Array a )
{
  return a->MAXSIZE;
}

void* ArrayGetData( Array a )
{
  return a->dataArray;
}

#undef __FUNCT__
#define __FUNCT__ "ArrayZero"
PetscErrorCode ArrayZero( Array a )
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscMemzero(a->dataArray, a->ELEMSIZE*a->MAXSIZE); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ArrayDelete1"
PetscErrorCode ArrayDelete1( Array a, int idx )
{
  /* Overwrites the deleted element of the array
   * with the last element: O(1) delete op
   */
  void* dest = a->dataArray + a->ELEMSIZE * idx;
  void* src  = a->dataArray + a->ELEMSIZE * (a->len-1);
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscMemcpy(dest,src,a->ELEMSIZE); CHKERRQ(ierr);
  a->len--;
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ArrayMap"
PetscErrorCode ArrayMap( Array a )
{
  int i;
//  PetscErrorCode ierr;

  PetscFunctionBegin;
  for ( i = 0; i < a->len; ++i) {
//    func( a->dataArray + a->ELEMSIZE * i );
  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ArrayCopy"
PetscErrorCode ArrayCopy( Array src, Array copy )
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = ArraySetSize(copy, ArrayLength(src) ); CHKERRQ(ierr);
  ierr = PetscMemcpy(copy->dataArray,src->dataArray, src->len*src->ELEMSIZE ); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
