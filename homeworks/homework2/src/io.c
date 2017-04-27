#include "struct.h"
#include <stdlib.h>
#include <stdio.h>

void print_physics_grid(physics_grid *P){
  int i;
  FILE *fp;
  fp = fopen("data.dat", "w");
  if (fp == NULL) {
    fprintf(stderr, "Can't open input file in.list!\n");
    exit(1);
  }
  for (i = 0; i < P->N; i++) {
    fprintf(fp,"%f %f %f %f %f \n", P->pos[i],P->rho[i],P->vel[i],P->pres[i],P->ene[i]);
  }
  fclose(fp);
}
