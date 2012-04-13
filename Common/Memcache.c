#include "Common.h"

struct _MemCache {
  size_t elemsize;
  size_t chunksize;
  Array cache;  // pointer array of free mem
  Array chunks; // pointer array of allocated chunks
};


#undef __FUNCT__
#define __FUNCT__ "MemCacheCreate"
PetscErrorCode MemCacheCreate( const char name[], size_t elemsize, size_t chunksize, MemCache *mc )
{
  MemCache mem;
  char tmp[256];
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscNew(struct _MemCache,&mem); CHKERRQ(ierr);
  sprintf(tmp, "%s_%s", name, "memcache");
  ierr = ArrayCreate( tmp, sizeof(void*), &mem->cache); CHKERRQ(ierr);
  sprintf(tmp, "%s_%s", name, "memchunk");
  ierr = ArrayCreate( tmp, sizeof(void*), &mem->chunks); CHKERRQ(ierr);
  mem->elemsize = elemsize;
  mem->chunksize = chunksize;
  *mc = mem;
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MemCacheDestroy"
PetscErrorCode MemCacheDestroy( MemCache mc )
{
  int i;
  void **chunks = ArrayGetData(mc->chunks);
  PetscErrorCode ierr;

  PetscFunctionBegin;
  for ( i = 0; i < ArrayLength(mc->chunks); ++i) {
    ierr = PetscFree(chunks[i]); CHKERRQ(ierr);
  }
  ierr = ArrayDestroy(mc->chunks); CHKERRQ(ierr);
  ierr = ArrayDestroy(mc->cache); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MemCacheAlloc"
PetscErrorCode MemCacheAlloc( MemCache mc, void *ptr )
{
  size_t elemsize = mc->elemsize;
  size_t chunksize = mc->chunksize;
  int len;
  int i;
  void **cache;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  len = ArrayLength(mc->cache);
  // if no more free mem to use
  if( len == 0 ) {
    void **chunk=0;
    //Get a new chunk ptr from chunkarray
    ierr = ArrayAppend(mc->chunks,&chunk); CHKERRQ(ierr);
    //Allocate chunk of new mem
    ierr = PetscMalloc(elemsize*chunksize,chunk); CHKERRQ(ierr);
    ierr = PetscMemzero(*chunk,elemsize*chunksize); CHKERRQ(ierr);
    //Add all free ptr to cache array
    ierr = ArraySetSize(mc->cache,chunksize); CHKERRQ(ierr);
    cache = ArrayGetData(mc->cache);
    //Mem added in reverse order since popped off last first
    for ( i = 0; i < chunksize; ++i) {
      cache[i] = *chunk + (chunksize-1-i)*elemsize;
    }
    const int s = ArrayLength(mc->chunks) * chunksize * elemsize / (1024*1024);
    ierr = PetscInfo1(0, "MemCache allocation: %d MB\n", s ); CHKERRQ(ierr);
  }
  //Pop off and return last mem location
  len = ArrayLength(mc->cache);
  ierr = ArrayGetP(mc->cache,len-1,ptr); CHKERRQ(ierr);
  ierr = ArraySetSize(mc->cache,len-1); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "MemCacheFree"
PetscErrorCode MemCacheFree( MemCache mc, void *ptr )
{
  void **pos=0;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscMemzero(ptr,mc->elemsize); CHKERRQ(ierr);
  ierr = ArrayAppend(mc->cache,&pos); CHKERRQ(ierr);
  *pos = ptr;
  PetscFunctionReturn(0);
}
