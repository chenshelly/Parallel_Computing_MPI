
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "mpi.h"
// #define SEQUENTIAL 
#define HEAVY  1000
#define SIZE       40
#define RADIUS 10
#define ERROR_CODE 999
void createValuesType(MPI_Datatype* dataType);
double calcHeavy(int x);

enum tags { WORK, STOP };
enum ranks { ROOT };
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


void workerProcess(MPI_Datatype valuesType)
{
   
    int chunk, tag;
    int index;
    MPI_Status status;
	double value=0;
    do
    {
        MPI_Recv(&chunk, 1, MPI_INT, ROOT, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        MPI_Recv(&index, 1, MPI_INT, ROOT, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        tag = status.MPI_TAG;
        for (size_t i = 0; i < chunk; i++)
           value += calcHeavy(index+i);

        MPI_Send(&value, 1, MPI_DOUBLE, ROOT, tag, MPI_COMM_WORLD);
    } while (tag != STOP);
}

void masterProcess(int numProcs, int chunk, MPI_Datatype valuesType,int *jobsDone)
{
    double time = MPI_Wtime();
    double sum=0;
    int remainder;
    int tag;
    double val=0;
    MPI_Status status;
for (int jobSent = 0; tag != STOP && jobSent < SIZE; jobSent += (numProcs - 1) * chunk)
{
    if (SIZE <= jobSent + (numProcs - 1) * chunk)
        {
            tag = STOP;
            chunk = (SIZE - jobSent) / (numProcs - 1);
            remainder = (SIZE - jobSent) % (numProcs - 1);
        }
    else tag=WORK;
    for (int workerId = 1; tag!=STOP&&workerId < numProcs; workerId++)
    {    
        printf("   %d-",*jobsDone);
        MPI_Send(&chunk,1,MPI_INT,workerId,tag,MPI_COMM_WORLD);
        MPI_Send(jobsDone,1,MPI_INT,workerId,tag,MPI_COMM_WORLD);
        MPI_Recv(&val,1,MPI_DOUBLE,workerId,tag,MPI_COMM_WORLD,&status);    
        sum+=val; 
        printf("%lf\n",val);
        *jobsDone+=chunk;
    }
}
	    if (remainder!=0)
    {
        printf("performing heavy for remainder \n");  //remainder part
        for (size_t i = SIZE-remainder; i < SIZE; i++)
            sum+=calcHeavy(i);     
    }
        printf("solution is :%lf\n",sum);       //print solution and runtime 
        printf("Run time: %lf\n", MPI_Wtime() - time);    
}

int main(int argc, char *argv[]) {
    int startIndex=0;
    int jobsDone=0;
    int myRank, numProcs, chunk = argc >= 2 ? atoi(argv[1]) : 2;
    printf("%d",numProcs);
    chunk=1;
    MPI_Init(&argc, &argv);
   
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs); 
    MPI_Datatype valuesType;
    createValuesType(&valuesType);

    //SEQUENTIAL PROGRAM
#ifdef SEQUENTIAL
    double time = MPI_Wtime();
    double sum=0;
    for (size_t i = 0; i < SIZE; i++)
    {
        for (size_t j = 0; j < SIZE; j++)
        {
            sum+=heavy(i,j);
        }
        
    }
    printf("solution is :%lf\n",sum);
    printf("Run time: %lf\n", MPI_Wtime() - time);
#endif    

    if (numProcs < 2)
    {
        MPI_Abort(MPI_COMM_WORLD, ERROR_CODE);
        exit(ERROR_CODE);
    }
  
    if (myRank == ROOT)
        masterProcess(numProcs, chunk, valuesType,&jobsDone);
    else
        workerProcess(valuesType);
   
    MPI_Finalize();
	return 0;
}

void createValuesType(MPI_Datatype* dataType)
{
    MPI_Type_contiguous(2, MPI_INT, dataType);
    MPI_Type_commit(dataType);
}

double calcHeavy(int x)
{
double sum=0;
for (size_t i = 0; i < SIZE; i++)
{
   sum+=heavy(x,i);
}

}
