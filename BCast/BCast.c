#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "mpi.h"

struct LossParams {
	double L;   // Loss value
	double *Q;  // Parameter list
};

int main(int argc, char* argv[]) {
	int rank;          /* rank of process      */
	int p;             /* number of processes  */
	int root = 0;      /* root proc            */
	MPI_Comm comm = MPI_COMM_WORLD;
	MPI_Status status; /* return status for    */
	int i;             /* for loop temp var    */
	int LQlen;         // number of parameters + 1
	size_t LQmem;      // size in bytes of LQ
	double* LQ;        // local [loss, p1, p2, ... ]
	double* LQarray=NULL;   // global [LQ1, LQ2, ... ]

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &p);

	//Read in param len;
	LQlen = 5;
	LQmem = LQlen * sizeof(double);
	LQ= (double*)malloc( LQmem );

	if( rank == root )
	{
		LQarray = (double*)malloc( LQmem * p );
		memset(LQarray,0,LQmem * p);
	}

	double start = MPI_Wtime();
	int j;
	for( j = 0; j < 1000; j++ )
	{
		// Run OPT
		srand(-rank);
		LQ[0] = rand() / 1000000.;
		int k;
		for( k = 1; k < LQlen; k++ )
		{
			LQ[k] = rank;
		}
		// Read LQ
		MPI_Gather(LQ,LQlen,MPI_DOUBLE,LQarray,LQlen,MPI_DOUBLE,root,comm);

		if( rank == root )
		{
			for( i = 0; i < p*LQlen; i+=LQlen )
			{
				if( LQarray[i] < LQ[0] )
				{
					memcpy(LQ,&LQarray[i],LQmem);
				}
			}
		}

		MPI_Bcast(LQ,LQlen,MPI_DOUBLE,root,comm);
		// Write LQ
		// run opt
	}
	if( rank == root )
	{
		double end = MPI_Wtime();
		printf("TIME: %f\n\n", end-start);
	}
    // Kill opt procs
	free( LQ );
	free( LQarray );
	MPI_Finalize();
	return 0;
} /* main */
