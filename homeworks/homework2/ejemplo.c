#include <stdio.h>
#include <stdlib.h>
#include <math.h>
int N = 10;
void init_to_zero(double * p, int n_points);
void mall(double ** matt);

int main(int argc, char const *argv[]) {
  double * mat;
  int i;

  mall(&mat);
  for (i = 0; i < N; i++) {
  //  printf(" fuera : %f \n", mat[i] );
  }
  return 0;
}

void mall(double ** matt) {
  int i;
  *matt=malloc(N * sizeof(double));
  printf("%lu\n",N * sizeof(double) );

//  init_to_zero(*matt,N);
  for (i = 0; i < 3; i++) {
   printf(" dentro : %f \n", *matt[i] );
  }
}
void init_to_zero(double * p, int n_points){
  int i;
  for(i=0;i<n_points;i++){
    p[i] = 0.0;
  }
}
