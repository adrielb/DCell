#include <stdlib.h>
#include <math.h>
#include "mpi.h"
#include "DCell.h"
#include "DWorld.h"

#define BSIZE 3000

void GetNeiRanks(MPI_Comm comm, int dims[], int coors[], int *neiRanks, int *numNei );
void LocalCellArrayAdd( DWorld world, DCell cell );
void LocalCellArrayRemove( DWorld world );
void SendDCell(DWorld world, int idx, int dst);
void SendEOF(DWorld world, int dst);

struct _DWorld {
	MPI_Comm comm; // MPI Cartesian comm
	double xs, ys; // Local start of domain (xs, ys)
	double gxs,gys;// Global start of domain (gxs,gys)
	double dx, dy; // Length of local domain

	int bufSize;  // MPI_Pack buffer size for DCell
	char *rPack;  // Receive buffer
	char *sPack;  // Send buffer

	// Crappy list data structure since glib not installed in cluster
	int localSize;         //Size of localCell array
	DCell localcells[BSIZE]; //Pointer array of cells
	int sendIdx[BSIZE];      //Index of localcells that will be sent
	int sendRank[BSIZE];     //Proc rank cell will be sent too
	int totSend;           //Total # of cells that will be sent

	int numNei;       //Number of non-null neighbor procs
	int neiRanks[9]; //Rank of all neighbors
};

void DWorldCreate(double gxs, double gxe, double gys, double gye, DWorld *world) {
	int i;
	struct _DWorld *w;
	w = malloc(sizeof(struct _DWorld));

	w->localSize = 0;
	for (i = 0; i < BSIZE; ++i) {
		w->localcells[i] = NULL;
	}

	int ndims = 2;
	int dims[3] = {0,0,0};
	int wrap[3] = {0,0,0};
	int reorder = 0;
	int size, rank;
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	MPI_Comm_size(MPI_COMM_WORLD,&size);
	MPI_Dims_create(size,ndims,dims);
	MPI_Cart_create(MPI_COMM_WORLD,ndims,dims,wrap,reorder,&w->comm);
	MPI_Comm_set_name(w->comm,"My Cart Comm");

	int coor[2];
	MPI_Cart_coords(w->comm,rank,ndims,coor);
	GetNeiRanks(w->comm, dims, coor, w->neiRanks, &w->numNei ); // in petsc: DAGetNeighbors()

	w->gxs = gxs;
	w->gys = gys;
	w->dx = (gxe - gxs) / dims[0];
	w->dy = (gye - gys) / dims[1];
	w->xs = gxs + coor[0] * w->dx;
//	w->xe = gxs + (coor[0]+1) * w->dx;
	w->ys = gys + coor[1] * w->dy;
//	w->ye = gys + (coor[1]+1) * w->dy;

	if( rank == 0 ) {
		printf("Dims: [%d, %d] at (%f, %f)\n", dims[0], dims[1], w->dx, w->dy);
	}

	//TODO: Make buffer size parameter adjustable
	w->bufSize = 50000;
	w->rPack = malloc(w->bufSize);
	w->sPack = malloc(w->bufSize);

	*world = w;
}

void DWorldDestroy( DWorld world ){
	MPI_Comm_free(&world->comm);
	free(world->rPack);
	free(world->sPack);
	free(world);
}

// Initialize to an equal density of cells per comm size
void DWorldInit( DWorld world ) {
	int i;
	DCell cell;
	int size;
	MPI_Comm_size(world->comm, &size);
	for (i = 0; i < (100/size); ++i) {
		DCellCreate(&cell);
		DCellInit(cell, world->xs,world->dx, world->ys,world->dy);
		LocalCellArrayAdd(world, cell);
	}
}

void VortexInABox(double x, double y, double *u, double *v);
void DWorldTimeStep(DWorld world, double dt )
{
	double x, y; // Coor of cell
	double u, v; // Velocity at coor
	DCell cell;
	int i;
	for (i = 0; i < world->localSize; ++i) {
		cell = world->localcells[i];
		DCellGetCoor(cell, &x, &y);
		VortexInABox(x,y,&u,&v);
		x = x + dt * u;
		y = y + dt * v;
		DCellSetCoor(cell, x,y);
	}
}

void VortexInABox(double x, double y, double *u, double *v) {
	const double PI = 3.141592653589793;
	*u = 2*cos(PI*y)*sin(PI*x)*sin(PI*x)*sin(PI*y);
	*v = -2*cos(PI*x)*sin(PI*x)*sin(PI*y)*sin(PI*y);
}

void DWorldBalance(DWorld world) {
	double x,y;
	DCell cell;
	int i, idx;
	int rank, myrank;
	int coor[2];
	MPI_Comm_rank(world->comm,&myrank);

	/* Calculate the rank the cell belongs too,
	 * check if it is equal to this proc's rank,
	 * if not, mark it in sendIdx and sendRank
	 */

	//Reset count to zero
	world->totSend = 0;

	for (i = 0; i < world->localSize; ++i) {
		cell = world->localcells[i];
		DCellGetCoor(cell, &x, &y);

		// Convert (x,y) to proc (i,j)
		coor[0] = floor( (x - world->gxs) / world->dx );
		coor[1] = floor( (y - world->gys) / world->dy );

		// Convert proc (i,j) to rank
		MPI_Cart_rank(world->comm,coor,&rank);

		// If the cell is off-proc,
		// record its index and its dst
		if( rank != myrank )
		{
		  // Just pack and send the cell, no need to index this
			idx = world->totSend;
			world->sendIdx[idx]  = i;
			world->sendRank[idx] = rank;
			world->totSend++;
		}
	}

	/* Perform the actual communication of the
	 * DCell objects. Interleave receives with
	 * sends.
	 */
	char *rPack = world->rPack;
	int bufSize = world->bufSize;
	int tag = 0;
	MPI_Request req;
	MPI_Status status;
	int countEOF = 0;
	int isEOF;

	i = 0;
	MPI_Barrier(world->comm);
	while( i <= world->totSend || countEOF < world->numNei )
	{
		if( countEOF < world->numNei  ) {
			MPI_Irecv(rPack,bufSize,MPI_PACKED,MPI_ANY_SOURCE,tag,world->comm,&req);
		}

		if( i < world->totSend ) {
			SendDCell(world, world->sendIdx[i], world->sendRank[i] );
//			LocalCellArrayRemove( world, world->sendIdx[i] );
			i++;
		}

		/* We are done sending, so alert
		 * all our neighbors with an EOF flag that this
		 * proc is complete
		 */
		if (i == world->totSend) {
			int j;
			for (j = 0; j < world->numNei; ++j) {
				SendEOF(world, world->neiRanks[j] );
			}
			i++;
		}

		if( countEOF < world->numNei ) {
			MPI_Wait(&req,&status);
			int pos = 0;
			MPI_Unpack(rPack,bufSize,&pos,&isEOF,1,MPI_INT,world->comm);
			if( isEOF == 1 ) {
				countEOF++;
			} else {
				DCell cell;
				DCellCreate(&cell);
				DCellUnpack(cell, rPack,bufSize,&pos,world->comm);
				LocalCellArrayAdd( world, cell );
			}
		}
	}

	/* NOW delete the cells that are sent
	 */
	LocalCellArrayRemove( world );
}

void DWorldSave(DWorld world, int time) {
	MPI_Status status;
	int i;
	int tag = 0;
	int root = 0; // Root proc that will write to disk
	int rank;
	int size;
	MPI_Comm_rank(world->comm, &rank);
	MPI_Comm_size(world->comm, &size);
//TODO: see if using proper MPI tags can eliminate this barrier
	MPI_Barrier(world->comm);
	if( rank == root ) {
		// Open file cells.time = [cell_1, cell_2, ...]
		FILE *fp;
		char filename[128];
//TODO: make directory location adjustable
		sprintf(filename,"data/cells.%d",time);
		fp = fopen(filename,"wb");
		// Root proc just saves cells to file locally
		int pos; // Position in packed buffer
		for (i = 0; i < world->localSize; ++i) {
			if( world->localcells[i] != NULL ) {
				pos = 0;
				DCellPack(world->localcells[i], world->rPack,world->bufSize,&pos,world->comm);
				DCellSave(world->rPack, fp );
			}
		}
		// Gather remote cells and write to file
		int countEOF = 0; // number of EOFs received
		int isEOF; // if received msg is an EOF when unpacked
		while( countEOF < size - 1 ) {
			MPI_Recv(world->rPack,world->bufSize,MPI_PACKED,MPI_ANY_SOURCE,tag,world->comm,&status);
			pos = 0;
			MPI_Unpack(world->rPack,world->bufSize,&pos,&isEOF,1,MPI_INT,world->comm);
			if( isEOF == 1 ) {
				countEOF++;
			} else {
				DCellSave( world->rPack + pos, fp );
			}
		}
		// Closes the file pointer
		fclose( fp );
	} else {
		// In a remote process, send all DCells to root.
		for (i = 0; i < world->localSize; ++i) {
			if( world->localcells[i] != NULL ) {
				SendDCell(world, i, root);
			}
		}
		// Signal root that no more DCells will come from this rank.
		SendEOF(world, root);
	}
}

void DWorldLoad(DWorld world, char *file ) {
	int i;
	int root = 0;
	int rank;
	MPI_Comm_rank(world->comm,&rank);
	int size;
	MPI_Comm_size(world->comm,&size);
	if( rank == root ) {
		int fileEOF = 1;
		while( fileEOF != 0) {
			// Open file for reading
			// Read (x,y)
			// Determine destination rank
			// Send to dst
		}
		// Let everybody know loading is complete
		for (i = 1; i < size; ++i) {
			SendEOF(world, i);
		}
	} else {
		int isEOF;
		while( 1 ) {
			// Blocking recv
			// Check EOF
			if( isEOF == 1 ) break;
			// unpack DCell
			// add to local list
		}
	}
}

/* Send DCell located at array index idx to
 * destination rank dst. Append eof = false
 * to pack buffer first.
 */
void SendDCell(DWorld world, int idx, int dst) {
	int tag = 0;
	int eof = 0;  // EOF IS FALSE
	int pos = 0;
	MPI_Pack(&eof,1,MPI_INT,world->sPack,world->bufSize,&pos,world->comm);
	DCellPack(world->localcells[idx],world->sPack,world->bufSize,&pos, world->comm);
	MPI_Ssend(world->sPack,pos,MPI_PACKED,dst,tag,world->comm);
}

/* Send eof = true to destination rank dst.
 * Signals to terminate receive requests
 */
void SendEOF(DWorld world, int dst) {
	int tag = 0;
	int eof = 1; // EOF IS TRUE
	int pos = 0;
	MPI_Pack(&eof,1,MPI_INT,world->sPack,world->bufSize,&pos,world->comm);
	MPI_Send(world->sPack,pos,MPI_PACKED,dst,tag,world->comm);
}

/* Add cell to end of array,
 * increment array size by one.
 */
void LocalCellArrayAdd( DWorld world, DCell cell )
{
	world->localcells[world->localSize] = cell;
	world->localSize++;
}

/* Based on the sendIdx list,
 * delete the DCell from memory,
 * then compact the list to remove NULL entries.
 */
void LocalCellArrayRemove( DWorld world ) {
	if( world->totSend == 0 ) return; //Nothing changed
	int i;
	int idx;
	for (i = 0; i < world->totSend; ++i) {
		idx = world->sendIdx[i];
		DCellDestroy( world->localcells[idx] );
		world->localcells[idx] = NULL;
	}
	idx = 0;
	DCell temp[BSIZE];
	for (i = 0; i < world->localSize; ++i) {
		if( world->localcells[i] != NULL ) {
			temp[idx] = world->localcells[i];
			idx++;
		}
	}
	world->localSize = world->localSize - world->totSend;
	for (i = 0; i < world->localSize; ++i) {
		world->localcells[i] = temp[i];
	}
}

/* Return the array of neighbor ranks, neiRanks,
 * and total number of non-null ranks, numNei
 */
void GetNeiRanks(MPI_Comm comm, int dims[], int coors[], int *neiRanks, int *numNeighbors ) {
	int c[2];
	int r, i, j;
	int numNei = 0;
	for (i = -1; i <= 1; i++) {
		for (j = -1; j <= 1; j++) {
			c[0] = coors[0] + i;
			c[1] = coors[1] + j;

			//Exclude yourself from neighbor list
			if( i==0 && j==0 )
				continue;
			/* If we are outside the process boundary,
			 * MPI_Cart_rank will exit with an error.
			 * So we have to check ourselves.
			 */
			if( c[0] == -1 || c[0] == dims[0] ||
				c[1] == -1 || c[1] == dims[1]  )
			{
				r = MPI_PROC_NULL;
			} else {
				MPI_Cart_rank(comm,c,&r);
				neiRanks[numNei] = r;
				numNei++;
			}
		}
	}
	*numNeighbors = numNei;
}
