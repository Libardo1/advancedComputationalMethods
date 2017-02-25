#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include <stdlib.h>

#define PI acos(-1.0)

// Parameters

double betta = 0.3;
int N = 1024;

// Define internal functions 

double *init(double *array);

int main(){
	
  	// MPI initialization 

  	MPI_Init(NULL,NULL);
  	int size, rank;
  	MPI_Comm_size(MPI_COMM_WORLD, &size);
  	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  	// We have to ckeck if the number of processors is allowed : # must be 2^n n = 1,2,3,4
 
	if(rank==0){
		if (size != 2 & size != 4 & size != 8 & size != 16 & size != 1){
     		printf("\n\n\n%d is not a number of processors allowed", size);
     		fprintf(stderr, "\nPlease run again with some of these numbers : 2, 4, 6, 8 \n\n");
     		abort();
     		MPI_Abort(MPI_COMM_WORLD, 1);
    	}
    	printf("\n\n\nWe have %d processors avaliables\n", size);
  	}
  
  	// Create the grids to store N points in the solid
  
  	double *grid = malloc(N*sizeof(double));
  	double *grid_new = malloc(N*sizeof(double));
  
  	// Initial conditions for the solid 
  
	grid = init(grid);
  	grid_new = init(grid);

	int T = 100*N;
  	MPI_Finalize();	
  	return(0);
}

double *init(double *array){ 
  
  //	int i;
//  	for(i=1; i<N-1; i++){
//    	array[i] = sin(2*PI*i/(double)(N-1));
//  }
	return array;
}
