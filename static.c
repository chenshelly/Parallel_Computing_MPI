
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "mpi.h"

#define HEAVY  1000
#define SIZE       40
#define RADIUS 10

// This function performs heavy computations, 
// its run time depends on x and y values
// DO NOT change the function
double heavy(int x, int y) {
	int i, loop;
	double sum = 0;

	if (sqrt((x - 0.75*SIZE)*(x - 0.75*SIZE) + (y - 0.25*SIZE)*(y - 0.25*SIZE)) < RADIUS)
		loop = 5*x*y;
	else
		loop = y + x;

	for (i = 0; i < loop*HEAVY; i++)
		sum += sin(exp(cos((double)i / HEAVY)));

	return  sum;
}

// Static code to be parallelized
int main(int argc, char *argv[]) {
	int X_start ,X_stop,Y_start=0 ,Y_stop=SIZE;
	int chunk;
	int rank;
    int numProcs;
	double local_sum = 0;
	double global_sum=0;
    double time;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

        if(numProcs < 2)
        {
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
            exit(EXIT_FAILURE);
            time = MPI_Wtime();
        }

    
	chunk=SIZE/numProcs; //number of tasks per processor
	X_start=chunk*rank; // processor's index start task
	X_stop = X_start+chunk;//processor's index stop task
	time = MPI_Wtime();
	while(X_start < X_stop)
	{
		for(int i=0 ; i < SIZE ; i++)
		{
			local_sum += heavy(X_start, i);
		}
        X_start++;
	}

	
	if(rank != 0)
    {
        MPI_Send(&local_sum, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    }
		
	else
	{
			double tempSum=0;
			global_sum=local_sum;
			for(int workerPrecs=1 ; workerPrecs < numProcs ; workerPrecs++)
			{
				MPI_Recv(&tempSum, 1, MPI_DOUBLE, workerPrecs, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				global_sum+=tempSum;
			}
		
			printf("The sum answer is:%e\n" , global_sum);
            printf("Run time: %lf\n", MPI_Wtime() - time);
	}

    MPI_Finalize();
    return 0;	
}
