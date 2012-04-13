#include "mpi.h"
#include "DCell.h"
#include <stdlib.h>

void CellBalancer(DWorld world, MPI_Comm comm) {
	int i;
	int nNei = 26;
	MPI_Status status[3*3*3];
	MPI_Request requests[3*3*3];
	int tag = MPI_ANY_TAG;
	int countSend[27];
	int countRecv[27];
	int totRecv;
	int totSend;

	/* Scatter how many cells this proc will send,
	 * and gather how many cells this proc will receive
	 */

	//Recieve counts from neighbors
	for( i = 0; i < nNei; i++ ) {
		MPI_Irecv(&countRecv[i],1,MPI_INT,MPI_ANY_SOURCE,tag,comm,&requests[i]);
	}

	//Scatter counts to neighbors
	for( i = 0; i < nNei; i++ ) {
		MPI_Isend(&countSend[i],1,MPI_INT,dest,tag,comm,&requests[i]);
	}
	MPI_Waitall(nNei,requests,status);

	for (i = 0; i < nNei; ++i) {
		totRecv += countRecv[i];
	}

	/* Perform the actual communication of the
	 * DCell objects.
	 */

	char *rPack = world->sPack;
	char *sPack = world->sPack;
	MPI_Request req;
	for( i = 0; i < max(totRecv, totSend); i++) {
		if( i < totRecv ) {
			MPI_Irecv(rPack,bufSize,MPI_PACKED,MPI_ANY_SOURCE,tag,comm,&req);
		}

		if( i < totSend ) {
			DCellPack(cell, buffer,bufSize, position, comm);
			MPI_Ssend(sPack,1,MPI_PACKED,dst,tag,comm);
			DCellDestroy(cell);
			Remove from cellArray
		}

		if( i < totRecv ) {
			MPI_Wait(req,status);
			DCell cell;
			DCellCreate(&cell);
			DCellUnpack(cell,rPack,bufSize,comm);
			Add to cell array
		}
	}

}
