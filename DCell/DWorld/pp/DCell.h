#ifndef DCELL_H_
#define DCELL_H_
#include "mpi.h"
#include <stdio.h>

typedef struct _DCell *DCell;

void DCellCreate( DCell *cell );
void DCellDestroy( DCell cell );
void DCellInit(DCell c, double xs, double dx, double ys, double dy);
void DCellPack(  DCell cell, void *buffer, int bufSize, int *pos, MPI_Comm comm );
void DCellUnpack(DCell cell, void *buffer, int bufSize, int *pos, MPI_Comm comm );
void DCellGetCoor(DCell cell, double *x, double *y);
void DCellSetCoor(DCell cell, double x, double y);
void DCellSave( void *buffer, FILE *fp);

#endif /* DCELL_H_ */
