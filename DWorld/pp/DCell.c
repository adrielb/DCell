#include "mpi.h"
#include "DCell.h"
#include <stdlib.h>
#include <math.h>

#define PI 3.14159265

struct _DCell {
	double x,y; // The (x,y) centroid of cell in global domain
	int ID;     // TODO: include cell ID in packed buffer
	int size;   // Array size of 2*n for cell outline
	double *coorX;  // [(x1,y1), (x2,y2), ...]
	double *coorY;
};

void DCellCreate( DCell *cell ) {
	DCell c;
	c = (struct _DCell*) malloc(sizeof(struct _DCell));
//TODO: create unique cell ID
//c->ID = 1000 * rank + ID;
	*cell = c;
}

void DCellDestroy( DCell cell ) {
	free( cell->coorX );
	free( cell->coorY );
	free( cell );
}

void DCellInit(DCell c, double xs, double dx, double ys, double dy ) {
	// randomly place cell in domain
	c->x = (rand() % 100 ) / 100.; // U[0,1]
	c->x = xs + c->x * dx; // [xs, xe] in local domain
	c->y = (rand() % 100 ) / 100.; // U[0,1]
	c->y = ys + c->y * dy; // [ys, ye] in local domain

	// randomly generate cells of memory size [1000-2000]
	int r = 100;
	c->size = rand() % r + r;
	c->coorX = (double *) malloc(c->size*sizeof(double));
	c->coorY = (double *) malloc(c->size*sizeof(double));

	double radius = 0.01 * (c->size / (double)r);
	int i;
	for( i = 0; i < c->size; i++ ) {
		c->coorX[i] = radius * cos(2*PI*i / c->size );
		c->coorY[i] = radius * sin(2*PI*i / c->size );
	}
}

// Pack DCell into buffer: [x, y, size, coorX, coorY]
void DCellPack(DCell cell, void *buffer, int bufSize, int *pos, MPI_Comm comm ) {
	MPI_Pack(&cell->x,1,MPI_DOUBLE,buffer,bufSize,pos,comm);
	MPI_Pack(&cell->y,1,MPI_DOUBLE,buffer,bufSize,pos,comm);
	MPI_Pack(&cell->size,1,MPI_INT,buffer,bufSize,pos,comm);
	MPI_Pack(cell->coorX,cell->size,MPI_DOUBLE,buffer,bufSize,pos,comm);
	MPI_Pack(cell->coorY,cell->size,MPI_DOUBLE,buffer,bufSize,pos,comm);
	//TODO: Record the size of 'pos' to figure out the typical size of the pack buffer
}

//Unpack [x,y,size,coorX,coorY] into DCell
void DCellUnpack(DCell cell, void *buffer, int bufSize, int *pos, MPI_Comm comm ) {
	MPI_Unpack(buffer,bufSize,pos,&cell->x,1,MPI_DOUBLE,comm);
	MPI_Unpack(buffer,bufSize,pos,&cell->y,1,MPI_DOUBLE,comm);
	MPI_Unpack(buffer,bufSize,pos,&cell->size,1,MPI_INT,comm);
	cell->coorX = (double *) malloc(cell->size*sizeof(double));
	cell->coorY = (double *) malloc(cell->size*sizeof(double));
	MPI_Unpack(buffer,bufSize,pos,cell->coorX,cell->size,MPI_DOUBLE,comm);
	MPI_Unpack(buffer,bufSize,pos,cell->coorY,cell->size,MPI_DOUBLE,comm);
}

void DCellGetCoor(DCell cell, double *x, double *y) {
	*x = cell->x;
	*y = cell->y;
}

void DCellSetCoor(DCell cell, double x, double y) {
	cell->x = x;
	cell->y = y;
}

/* reuse packed buffer binary layout for writing to disk
 */
void DCellSave( void *buffer, FILE *fp) {
	//Do some pointer math to extract the size
	int *size = buffer + 2*sizeof(double);
	// totalBytes = (x, y, coorX, coorY) * 8 + size * 4
	int totalBytes = sizeof(double)*(2+2*(*size)) + sizeof(int);
	// write entire buffer to disk
	fwrite(buffer,sizeof(char),totalBytes,fp);
}
