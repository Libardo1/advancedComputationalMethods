#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>

// Some parameters to use in the code
int L = 5;
int d = 1, l = 2, V0 = 100, m = 8, N;
double h;

// Functions 
int transformer(int i, int j);
double *init(int x0, int x1, int y0, int y1, double *array);
void ind2sub(int *i, int *j, int a, int b, int ind);
double *init_alt(int row, int col, double *array);

int main(int argc, char** argv){

   	MPI_Init(NULL, NULL);	
   		
	// Get the number of processes  
   	int world_size;
   	MPI_Comm_size(MPI_COMM_WORLD, &world_size); 
   	
   	// Get the rank of the process
   	int world_rank;
   	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
   	
   	//We have to ckeck if the number of processors is allowed : # must be 2^n, n = 1,2,3,4
	if(world_rank==0)
	{
		if (world_size != 2 & world_size != 4 & world_size != 8 & world_size != 16 & world_size != 1)
		{
			printf("\n\n\n%d is not a number of processors allowed", world_size);
			fprintf(stderr, "\nPlease run again with some of these numbers : 2, 4, 6, 8 \n\n");
			abort();
			MPI_Abort(MPI_COMM_WORLD, 1);
		}
		printf("\n\n\nWe have %d processors avaliables\n", world_size);
	} 
	
//	int up, down, left, right, x0, x1, y0, y1, i, j, k, n=0;
//  double average;
//  N = 2*m*m;
//	h = L/(double)m;   		
   	
//	printf("m : %d, l %d, h %f, d %d \n",m,l,h,d);
	
//	x0 = m/2 - l/(h*2) - 1;
//	x1 = m/2 + l/(h*2) - 1;
//	y0 = m/2 - d/(h*2) - 1;
//	y1 = m/2 + d/(h*2) - 1;
	
//	printf("x0 %d, x1 %d, y0 %d, y1 %d\n",x0,x1,y0,y1);
	
	// Create the total array in the root process (0)
//	double V, V_new;
//	if(world_rank == 0)
//	{	
//		double *V 		= malloc(m*m*sizeof(double));
//		double *V_new 	= malloc(m*m*sizeof(double));

//		V = init_alt(m, m, V);

//		V 		= init(x0, x1, y0, y1, V);		
//		V_new 	= init(x0, x1, y0, y1, V);
//	}
	
	
	// For each process create a sub-array to hold a part 
	// of the total array in each iteration
//	int sub_m = (m/world_size);
//	int sub_m_c 	= (sub_m + 2) * m;	// Size of each subgrid for n = 1,...,n-1
//	int sub_m_b 	= (sub_m + 1) * m; 	// Size of each subgrid for n = 0,n

//	double *sub_V, *sub_V_new;
		
//	if(world_rank == 0 | world_rank == world_size-1)
//	{
//		double *sub_V 		= malloc(sub_m_b*sizeof(double));		
//		double *sub_V_new 	= malloc(sub_m_b*sizeof(double));
//	}
//	else
//	{
//		double *sub_V 		= malloc(sub_m_c*sizeof(double));		
//		double *sub_V_new 	= malloc(sub_m_c*sizeof(double));
//	}	

//	MPI_Barrier( MPI_COMM_WORLD );

	// Set the initial and final limits to cut data in each process
	// Values are linear index oriented in column major order index E (1,m*m)
//	int initial, final;
//	if (world_rank == 0)
//	{
//		initial = 1;
//		final 	= (sub_m + 1) * m -1;
//	}
//	else if (world_rank == world_size-1 )
//	{
//		initial = (sub_m * world_rank - 1) * m  ;
//		final 	= m*m -1;	
//	}		
//	else
//	{
//		initial = (sub_m * world_rank - 1) * m ;	
//		final 	= (sub_m * (world_rank +1) + 1) * m - 1 ; 
//	}
//	int size 	= final - initial + 1;

//	if (1 == 1)			// TEST
//	{
//		printf("I am Process %d with initial: %d and final: %d and size :%d\n",world_rank, initial, final, size);
//	}
//	if (size != sizeof(sub_V))
//	{
//		printf("ERROR: bad shape arrays :: size : %d vs sub_V : %d\n ", size, (int)sizeof(sub_V));
//	}
//
//	if (world_rank == 0)
//	{	
//		for (int pro = 0; pro<world_size; pro++)
//		{				
//			printf("Sending data to process %d\n", pro);
//			if(pro == 0){
//				double array[(sub_m +1) * m] = V[0:(sub_m + 1) * m];
//				MPI_Send(&V[1], 1, MPI_DOUBLE, pro, 0, MPI_COMM_WORLD);	
//			}
//			else if(pro == world_size-1){
//				double array[m*m -1] = V[0:(sub_m + 1) * m];
//				MPI_Send(&V[1], 1, MPI_DOUBLE, pro, 0, MPI_COMM_WO
//			}
//			else{
//				double array[(sub_m +1) * m] = V[0:(sub_m + 1) * m];
//				MPI_Send(&V[1], 1, MPI_DOUBLE, pro, 0, MPI_COMM_WO
//			}
//		}
//	}
//	MPI_Recv(&sub_V, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//	printf("I am Process %d with initial: %d and final: %d \n",world_rank, initial, final);
	
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

double *init_alt(int row, int col, double *array){	
	
	int i,j;
	for(j = 0; j < col; j++){
		for(i = 0; i < row; i++){
//			printf("linear: %d, i: %d, j: %d\n", transformer(i, j), i, j);
			array[transformer(i, j)] = transformer(i, j);
		}
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