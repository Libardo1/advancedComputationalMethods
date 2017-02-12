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
   	
   	//We have to ckeck if the number of processors is allowed : # must be 2^n n = 1,2,3,4
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
	
	int up, down, left, right, x0, x1, y0, y1, i, j, k, n=0;
   	double average;
   	N = 2*m*m;
   	h = L/(double)m;   		

//	printf("m : %d, l %d, h %f, d %d \n",m,l,h,d);

	x0 = m/2 - l/(h*2) - 1;
	x1 = m/2 + l/(h*2) - 1;
	y0 = m/2 - d/(h*2) - 1;
	y1 = m/2 + d/(h*2) - 1;
	
//	printf("x0 %d, x1 %d, y0 %d, y1 %d\n",x0,x1,y0,y1);

// Create the total array in each process
	double *V 		= malloc(m*m*sizeof(double));
	double *V_new 	= malloc(m*m*sizeof(double));
	V 		= init(x0, x1, y0, y1, V);
	V_new 	= init(x0, x1, y0, y1, V);
	
	double *l_send 	= malloc(m*sizeof(double));
	double *r_send 	= malloc(m*sizeof(double));
	
	double *l_recv 	= malloc(m*sizeof(double));
	double *r_recv 	= malloc(m*sizeof(double));
	
//	if (world_rank == 0)
//	{ int h,k;
//		for(i = 0; i < 4; i++){
//			for(j = 0; j < 4; j++){
//				ind2sub(&h, &k, 4, 4,transformer(i,j));
//				printf("i : %d, j : %d, linear : %d, i_p : %d, j_p : %d\n", i,j,transformer(i,j),h,k);
//			}
//		}
//	}

	MPI_Barrier( MPI_COMM_WORLD );
	int sub_m = (m/world_size);				// Size of each sub part of the total grid in each process

	int initial, final;
	if (world_rank == 0)
	{
		initial = 1;
		final 	= (sub_m + 1);
	}
	else if (world_rank == world_size-1)
	{
		initial = (sub_m * world_rank - 1)   ;
		final 	= m -1;	
	}		
	else
	{
		initial = (sub_m * world_rank - 1)  ;	
		final 	= (sub_m * (world_rank +1) + 1) ; 
	}
	
	int size 	= final - initial + 1;
	
			
	while (n < N)
	{	
		// We divide the grid in equal-size parts according to the number of processes 	
		k = 0;
		for(j = initial; j <= final; j++)
		{
			for(i = 1; i < m; i++)
			{				
				if( k < sub_m)
				{
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
				else{
					i = m;
					j = final+1;
				}
				k++;
			}
		}
			
		MPI_Barrier( MPI_COMM_WORLD );
		
		k = 0;
		for(j = initial; j <= final; j++)
		{
			for(i = 1; i < m; i++)
			{				
				if( k < sub_m)
				{		
					V[transformer(i,j)] = V_new[transformer(i,j)];	
				}
				else{
					i = m;
					j = final+1;
				}
				k++;
			}
		}
				
				
		// Preparing data to send between process
		if(world_rank == 0)
		{
			for (i = 0; i < m; i++)
			{
				r_send[i] = V_new[transformer(i,sub_m-1)];
			}
			MPI_Send(&r_send, 1, MPI_FLOAT, world_rank+1, 1, MPI_COMM_WORLD);
		}
		else if(world_rank == world_size-1)
		{
			for (i = 0; i < m; i++)
			{
				l_send[i] = V_new[transformer(i,sub_m*world_rank+1)];
			}	
			MPI_Send(&l_send, 1, MPI_FLOAT, world_rank-1, 0, MPI_COMM_WORLD);

		}
		else
		{	
			for (i = 0; i < m; i++)
			{
				l_send[i] = V_new[transformer(i,sub_m*world_rank+1)];
				r_send[i] = V_new[transformer(i,sub_m*(world_rank+1)-1)];
			}	
			MPI_Send(&r_send, 1, MPI_FLOAT, world_rank+1, 1, MPI_COMM_WORLD);
			MPI_Send(&l_send, 1, MPI_FLOAT, world_rank-1, 0, MPI_COMM_WORLD);
		}
		
		// Recieve information between process
		if(world_rank == 0)
		{
			MPI_Recv(&r_recv, 1, MPI_FLOAT, world_rank+1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			for (i = 0; i < m; i++)
			{
				V[transformer(i,sub_m-1)] = r_recv[i];
			}
		}
		else if(world_rank == world_size-1)
		{
			MPI_Recv(&l_recv, 1, MPI_FLOAT, world_rank-1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			for (i = 0; i < m; i++)
			{
				V[transformer(i,sub_m*world_rank+1)] = l_recv[i];
			}
		}
		else
		{	
			MPI_Recv(&l_recv, 1, MPI_FLOAT, world_rank-1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			MPI_Recv(&r_recv, 1, MPI_FLOAT, world_rank+1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			for (i = 0; i < m; i++)
			{
				V[transformer(i,sub_m*world_rank+1)] 	 = l_recv[i];
				V[transformer(i,sub_m*(world_rank+1)-1)] = r_recv[i];
			}
		}
		
		n += 1;
		if( world_rank == 0)
		{
			printf("n : %d\n",n);
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
  	
  	*i = k;
  	*j = h;
}