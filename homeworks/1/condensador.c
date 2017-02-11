#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>

//int L = 5;
int d = 1, l = 2, V0 = 100, m = 128, N;
double h = 0.02;

int transformer(int i, int j);
double *init(int x0, int x1, int y0, int y1, double *array);
void ind2sub(int *i, int *j, int a, int b, int ind);

int main(int argc, char** argv){

   	MPI_Init(NULL, NULL);
	
	// First we check the world in we are
	
	// Get the number of processes  
   	int world_size;
   	MPI_Comm_size(MPI_COMM_WORLD, &world_size); 
   	
   	// Get the rank of the process
   	int world_rank;
   	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
   	
   	//We have to ckeck if the number of processors is allowed : # must be 2^n n = 1,2,3,4
	if(world_rank==0)
	{
		if (world_size != 2 & world_size != 4 & world_size != 8 & world_size != 16 & world_size != 1)
		{
			printf("\n %d is not a number of processors allowed", world_size);
			fprintf(stderr, "\nPlease run again with some of these numbers : 2, 4, 6, 8 \n\n");
			abort();
			MPI_Abort(MPI_COMM_WORLD, 1);
		}
		printf("We have %d processors avaliables\n", world_size);
	} 
	
	int up, down, left, right, x0, x1, y0, y1, i, j, k, n=0;
   	double average;

	// Metodo de relajación 
	
//	m = L/h; 
	
	N = 2*m*m;
	
	x0 = m/2 - l/(h*2) - 1;
	x1 = m/2 + l/(h*2) - 1;
	y0 = m/2 - d/(h*2) - 1;
	y1 = m/2 + d/(h*2) - 1;
	
	double *V = malloc(m*m*sizeof(double));
	double *V_new = malloc(m*m*sizeof(double));

	V = init(x0, x1, y0, y1, V);
	V_new = init(x0, x1, y0, y1, V);
		
	// Size of the sub-grids

	int m_k;
	m_k = (m*m)/world_size;
	if(world_rank==0)
	{
		printf("\nThe size of each sub grid is %d from a total of %d. \n",m_k, m*m);
	}
	
//	MPI_Barrier( MPI_COMM_WORLD );
	
	if (world_rank == 0 )
	{
//			ind2sub(&i, &j, 2, 2, 4);
//			printf("i %d, j %d\n", i, j );
	}
	
	if (world_rank == 1)
	{
//		V[1]=1.0;
	}
	
  	while (n < N)
	{	
		// We divide the grid in equal-size parts according to the number of processes 	
		for(k = 1; k < m_k+1; k++)
		{
			ind2sub(&i, &j, m, m, world_rank * m_k + k );
			if (i != 0 & j!=0 & i!=m & j!=m){
				up 		= transformer(i-1, j);
				down 	= transformer(i+1, j);
				left 	= transformer(i, j-1);
				right 	= transformer(i, j+1);
				if (!(j >= x0 && j <= x1 && i == y0) && !(j >= x0 && j <= x1 && i == y1))
				{	
					average = (V[up] + V[down] + V[left] + V[right])/4;
					V_new[transformer(i,j)] = average;
				}
			}
		}
		MPI_Barrier( MPI_COMM_WORLD );
		
		for(k = 1; k < m_k+1; k++)
		{
			if (i != 0 & j!=0 & i!=m & j!=m)
			{
				V[world_rank * m_k + k] = V_new[world_rank * m_k + k];
			}
		}

		n += 1;

		if(world_rank==0)
		{
			printf("n		   %d from process %d \n", n, world_rank);
		}
	}
	
	MPI_Finalize();	
	return(0);
}

int transformer(int i, int j){	
	
	return j*m + i;
}

double *init(int x0, int x1, int y0, int y1, double *array){	
	
	int a;
	for(a = x0; a <= x1; a++){
		array[transformer(y0, a)] = V0/2;
		array[transformer(y1, a)] = -V0/2;
	}
	return array;
}

void ind2sub(int *i, int *j, int a, int b, int ind){

	int h,k;
	 
  	h = (int)((ind/(double)a)-0.001);
  	h ++;
  	k = ind - (a*(h-1));
  	
  	*i = k-1;
  	*j = h-1;
}