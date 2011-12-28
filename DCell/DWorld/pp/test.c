#include "mpi.h"
#include "DWorld.h"
#include <stdio.h>

const int MAXITER = 400;
const double dt = 0.01;


int main(int argc, char* argv[]) {
	MPI_Init(&argc,&argv);

	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	double startBal, avgBal=0, startSave, avgSave=0;

	// Create a world that is [0,1]x[0,1]
	DWorld world;
	DWorldCreate(0,1,0,1, &world);

	// Randomly initialize cells in domain
	DWorldInit(world);

	// Write cells to file on root node
	DWorldSave(world, 0);

	// The main simulation loop
	int i;
	for (i = 1; i < MAXITER; ++i) {

		if( rank == 0 ) printf("i: %d\n", i);

		// Moves the cells around the domain
		DWorldTimeStep(world, dt);

		// Communicate off-proc cells
		startBal = MPI_Wtime();
		DWorldBalance(world);
		avgBal += ( MPI_Wtime() - startBal ) / MAXITER;

		// Save state for visualization
		startSave = MPI_Wtime();
		DWorldSave(world, i);
		avgSave += ( MPI_Wtime() - startSave ) / MAXITER;
	}

	if( rank == 0 ) {
		printf("Balancing Time: %f\n", avgBal);
		printf("Saving Time: %f\n", avgSave);
	}

	DWorldDestroy(world);
	MPI_Finalize();
	return 0;
}
