/* NOT USED
 * PREVIOUS IDEA OF COMMUNICATING RECV COUNT
 * PRIOR TO MPI PACKED COMMUNICATION
 * NEW CODE MERGED BOTH INTO MPI-PACKED
 */

void CountSends(DWorld world) {
	double x,y;
	DCell cell;
	int i, idx;
	int rank, myrank;
	int coor[2];
	MPI_Comm_rank(world->comm,&myrank);

	//Reset counts to zero
	world->totSend = 0;
	for (i = 0; i < 27; ++i) {
		world->countSend[i] = 0;
	}
	for (i = 0; i < world->localSize; ++i) {
		cell = world->localcells[i];
		DCellGetCoor(cell, &x, &y);

		// Convert (x,y) to proc (i,j)
		coor[0] = floor( (x - world->gxs) / world->dx );
		coor[1] = floor( (y - world->gys) / world->dy );

		// Convert proc (i,j) to rank
		MPI_Cart_rank(world->comm,coor,&rank);

		// If the cell is off-proc,
		// record its index, accumulate send count
		if( rank != myrank )
		{
			idx = world->totSend;
			world->sendIdx[idx]  = i;
			world->sendRank[idx] = rank;
			world->totSend++;
			idx = world->rankToSendCountIdx[rank];
			world->countSend[idx]++;
		}
	}
}

void DWorldBalance(DWorld world, MPI_Comm comm) {
	int i;
	int nNei = 26;
	MPI_Status status[3*3*3];
	MPI_Request requests[3*3*3];
	int tag = 0;
	int countSend[27];
	int countRecv[27];
	int totRecv;
	int totSend;

	CountSends(world);

	/* Scatter how many cells this proc will send,
	 * and gather how many cells this proc will receive
	 */

	//Recieve counts from neighbors
	for( i = 0; i < nNei; i++ ) {
		MPI_Irecv(&countRecv[i],1,MPI_INT,MPI_ANY_SOURCE,tag,comm,&requests[i]);
	}

	//Scatter counts to neighbors
	for( i = 0; i < nNei; i++ ) {
		MPI_Isend(&countSend[i],1,MPI_INT,world->neiRanks[i],tag,comm,&requests[i]);
	}
	MPI_Waitall(nNei,requests,status);

	for (i = 0; i < nNei; ++i) {
		totRecv += countRecv[i];
	}

	/* Perform the actual communication of the
	 * DCell objects. Interleave receives with
	 * sends.
	 */

	DCell cell;
	char *rPack = world->sPack;
	char *sPack = world->sPack;
	int bufSize = world->bufSize;
	MPI_Request req;
	int j=0; // indexes total send count
	int maxI = totRecv > totSend ? totRecv : totSend;
	for( i = 0; i < maxI; i++) {
		if( i < totRecv ) {
			MPI_Irecv(rPack,bufSize,MPI_PACKED,MPI_ANY_SOURCE,tag,comm,&req);
		}

		if( i < totSend ) {
			cell = world->localcells[world->sendIdx[j]];
			j++;
			int position;
			DCellPack(cell,sPack,bufSize,&position, comm);
			MPI_Ssend(sPack,position,MPI_PACKED,world->sendRank[j],tag,comm);
			LocalCellArrayRemove( world, world->sendIdx[j] );
		}

		if( i < totRecv ) {
			MPI_Wait(&req,status);
			int pos = 0;
			// Unpack isEOF
			exit(0);
			DCell cell;
			DCellCreate(&cell);
			DCellUnpack(cell, rPack,bufSize,&pos,comm);
			LocalCellArrayAdd( world, cell );
		}
	}

}
