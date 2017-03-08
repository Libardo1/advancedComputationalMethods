#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define PI acos(-1.0)

// Parameters

int N = 64;
double delta_t  = 5E-3;
double betta = 1.0;

// Define internal functions

double* init(double *array);
void F(double * F_array, double * x);

// Main program

int main(){

  	int T = 5*pow(N,2.2);
    int i,j,k,t,n,cont = 0;
    int P = 1000, K = 3;
    double Q, Qp, omega = 1.0;

  	// Create the grids to store N points in the solid

  	double* x = malloc(N*sizeof(double));
  	double* x_new = malloc(N*sizeof(double));
  	double* v = malloc(N*sizeof(double));
  	double* v_new = malloc(N*sizeof(double));
  	double* F_grid = malloc(N*sizeof(double));
  	double* F_grid_new = malloc(N*sizeof(double));

  	// Create the total grid to store data

    double **energy = (double **) malloc(K * sizeof(double *));
    for(i=0; i<K; i++){
		    energy[i] = (double *) malloc(P * sizeof(double *));
	  }

  	// Initial conditions for the solid
    // Position
	  x = init(x);
  	x_new = init(x_new);
    // Velocity & Acceleration
    for(n=0;n<N;n++){
      v[n] = 0.0;
      v_new[n] = 0.0;
      F_grid[n] = 0.0;
      F_grid_new[n] =0.0;
    }

    // **********************************************	//
	  // 							Leapfrog method										//
	  // **********************************************	//
	  // We use the Leapfrog method expressed in x, v & a quantities with integer steps
    for(t = 0; t < T; t++){
      F(F_grid,x);
      for(n = 1; n < N-1; n++){
        x_new[n] = x[n] + ( v[n] * delta_t ) + (0.5 * F_grid[n] * pow(delta_t,2.0) );
      }
      F(F_grid_new,x_new);
      for(n = 1; n < N-1; n++){
        v_new[n] = v[n] + 0.5 * (F_grid[n] + F_grid_new[n]) * delta_t;
      }

      // Update
      x = x_new;
      v = v_new;

//  		printf("iter: %d \n",t%N);

      if(t%(T/P) == 0){
        // Calculate the general position
        Q = 0.0; Qp = 0.0;
        for ( k = 1; k <= K; k++) {
          for ( n = 0; n < N; n++) {
            Q = Q + x[n]*sin((PI*k*n)/(double)N);
            Qp = Qp + v[n]*sin((PI*k*n)/(double)N);
          }
          Q = Q * sqrt(2.0/(double)N);
          Qp = Qp * sqrt(2.0/(double)N);
          energy[k-1][cont] = 0.5 * (pow(Q,2.0)+pow(omega*Qp,2.0));
        }
        cont++;
      }
    }

    // w2?
    // Print results
    for(i = 0; i<K;i++){
      for(j = 0; j<P;j++){
        printf("%f \n", energy[i][j]);
      }
    }
  	return(0);
}

// Auxiliar functions

double *init(double *array){
    int i;
  	array[0] = 0.0;
  	array[N-1] = 0.0;

	  for(i=1; i<N-1; i++){
		    array[i] = sin(PI*(double)i/(double)(N-1));
	  }
	  return array;
}

void F(double * F_array, double * x){
	int n;
	for(n = 1; n < N-1; n++){
		F_array[n] = (x[n+1] - 2.0 * x[n] + x[n-1]) + betta * (pow((x[n+1] - x[n]),3.0) - pow((x[n] - x[n-1]),3.0) );
	}
}
