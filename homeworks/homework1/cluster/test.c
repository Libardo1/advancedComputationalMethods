#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <omp.h>
#include <time.h>

int main() {
  int N = 1;
  double* vector = malloc(N*sizeof(double));
  omp_set_num_threads(3);
  #pragma omp parallel
  printf("Hello from thread %d, nthreads %d\n", omp_get_thread_num(), omp_get_num_threads());

  #pragma omp paralle for
    for (int i = 0;i<N;i++){
      for (int j = 0;i<N;j++){
          printf("%d from thread %d \n", 1, omp_get_thread_num());
          vector[i] = (double)1;
      }
    }

  return (0);
}
