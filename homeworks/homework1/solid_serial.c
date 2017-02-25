#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define PI acos(-1.0)

// Parameters

double betta = 0.3;
int N = 100;
double delta_t  = 10E-3;

// Define internal functions 

double* init(double *array);
void F(double * F_array, double * x);

int main(){
  
  	int T = 100*N;
	int i,j,k,t,n,cont = 0;
  
  	// Create the grids to store N points in the solid
  
  	double* x = malloc(N*sizeof(double));
  	double* x_new = malloc(N*sizeof(double));
  	double* v = malloc(N*sizeof(double));
  	double* v_new = malloc(N*sizeof(double));
  	double* F_grid = malloc(N*sizeof(double));
  	double* F_grid_new = malloc(N*sizeof(double));
  	
  	// Create the total grid to store data
  	
	double **DATA = (double **) malloc(N * sizeof(double *));
	for(i=0;i<N;i++){
		DATA[i] = (double *) malloc(N * sizeof(double *));
	}
  	
  	// Initial conditions for the solid 
  
	x = init(x);
  	x_new = init(x_new);
	
	for(n=0;n<N;n++){
		v[n] = 0.0;
		v_new[n] = 0.0;
		F_grid[n] = 0.0;
		F_grid_new[n] =0.0; 
	}
	
	// ***************************************************************************	//
	// 							Leapfrog method										//
	// ***************************************************************************	//
	// We use the Leapfrog method expressed in x, v & a quantities with integer steps 
	
	for(t = 0; t < T; t++){
		
		// Update
		
		x = x_new;
		v = v_new;
		
		F(F_grid,x);
		for(n = 1; n < N-1; n++){
			x_new[n] = x[n] + v[n] * delta_t + 0.5 * F_grid[n] * pow(delta_t,2);
		}	
		F(F_grid_new,x);
		for(n = 1; n < N-1; n++){
			v_new[n] = v[n] + 0.5 * (F_grid[n] + F_grid_new[n]) * delta_t;
		}
	
		if(t%N == 0){
			for(i = 0; i<N;i++){
				DATA[i][cont] = x_new[i];
			}
			cont++;
		}				
	}
	
	
	// print results
	
	for(i = 0; i<N;i++){
		for(j = 0; j<N;j++){
			printf("%f \n", DATA[i][j]);
		}	
	}
  	return(0);
}

double *init(double *array){ 
  	int i;
  	array[0] = 0.0;
  	array[N-1] = 0.0;
	for(i=1; i<N-1; i++){
		array[i] = sin(2*PI*i/(double)(N-1));
	}
	return array;
}

void F(double * F_array, double * x){
	int n;
	for(n = 1; n < N-1; n++){
		F_array[n] = (x[n+1] - 2.0 * x[n] + x[n-1]) + 
				betta * ( pow(x[n+1] - x[n],3.0) - pow(x[n] - x[n-1],3) );
	}
}