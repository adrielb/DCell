#ifndef ARRAY_H_
#define ARRAY_H_

#include "petsc.h"

typedef struct _Array *Array;

PetscErrorCode ArrayCreate( const char name[], int elemSize, size_t maxSize, Array *array );
PetscErrorCode ArrayDestroy( Array a );
PetscErrorCode ArraySetSize( Array a, int size );
PetscErrorCode ArrayAppend( Array a, void *elem );
PetscErrorCode ArrayAppendPtr( Array a, void *elem );
PetscErrorCode ArrayGet( Array a, int i, void *elem );
PetscErrorCode ArrayGetP( Array a, int i, void *elem );
PetscErrorCode ArrayWrite( Array a, const char *name, int t );
PetscErrorCode ArrayZero( Array a );
PetscErrorCode ArrayDelete1( Array a, int idx );
PetscErrorCode ArrayCopy( Array src, Array copy );
PetscErrorCode ArrayDuplicate( Array a, Array *newarray );

int   ArrayLength( Array a );
void* ArrayGetData( Array a );

#endif /* ARRAY_H_ */
