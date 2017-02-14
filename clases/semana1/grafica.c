#include <stdio.h>

int main(void){
  unsigned int id;
  int min;
  int max;
  float x;

  unsigned long int id2;
  long int min2;
  double x2;

  min = -12;
  max = 12;
  
  for (id=min;id<max+20;id++){
    x = (*(float*)&1E(float*)id);
    printf("%.25e\n", x);
  }

  return 0;
}
