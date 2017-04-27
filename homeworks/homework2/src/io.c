#include <stdlib.h>
#include <stdio.h>
#include "struct.h"
#include "calibration.h"
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

void print_calibration(){
  FILE *fp;
  fp = fopen("calibration.dat", "w");
  if (fp == NULL) {
    fprintf(stderr, "Can't open input file in.list!\n");
    exit(1);
  }
  fprintf(fp,"%f %f %f %f %f %f %f %f %f %f %f \n", T,CFL,X0,RL,PL,V0,RR,PR,DELTA_X,GAMMA,LEN_TUBE);
  fclose(fp);
}
