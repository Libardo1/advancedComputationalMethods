#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "calibration.h"
#include "struct.h"

double energy(double p,double r,double v){
  return (p/(GAMMA-1)) + (r*pow(v,2.0))/2.0;
}
double presion(double e, double r,double v){
    return (GAMMA-1)*(e - r *pow(v,2.0)/2.0);
}
void init_to_zero(FLOAT *p, int n_points){
  int i;
  for(i=0;i<n_points;i++){
    p[i] = 0.0;
  }
}
physics_grid * create_physics_grid(void){
  physics_grid *G;
  if(!(G = malloc(sizeof(physics_grid)))){
    fprintf(stderr, "Problem with data allocation\n");fflush(stdout);
    exit(0);
  }
  G->N=0;
  G->L=0.0;
  G->delta_x=0.0;
  G->pos=NULL;
  G->pres=NULL;
  G->rho=NULL;
  G->vel=NULL;
  G->ene=NULL;
  return G;
}
U_grid * create_U_grid(void){
  U_grid *G;
  if(!(G = malloc(sizeof(U_grid)))){
    fprintf(stderr, "Problem with data allocation\n");fflush(stdout);
    exit(0);
  }
  G->N=0;
  G->U_1=NULL;
  G->U_2=NULL;
  G->U_3=NULL;
  G->U_1_star=NULL;
  G->U_2_star=NULL;
  G->U_3_star=NULL;
  G->U_1_new=NULL;
  G->U_2_new=NULL;
  G->U_3_new=NULL;

  return G;
}
F_grid * create_F_grid(void){
  F_grid *G;
  if(!(G = malloc(sizeof(F_grid)))){
    fprintf(stderr, "Problem with data allocation\n");fflush(stdout);
    exit(0);
  }
  G->N=0;
  G->F_1=NULL;
  G->F_2=NULL;
  G->F_3=NULL;
  G->F_1_star=NULL;
  G->F_2_star=NULL;
  G->F_3_star=NULL;
  return G;
}
void init_problem(physics_grid *P, U_grid *U, F_grid *F){
  int i;
  P->L = LEN_TUBE;
  P->delta_x = DELTA_X;

  P->N = (int)(P->L/P->delta_x);
  U->N=P->N;
  F->N=P->N;

  if(!(P->pos=malloc(P->N * sizeof(FLOAT)))){
    fprintf(stderr, "Problem with pressure allocation");
    exit(1);
  }
  for ( i = 0; i < P->N; i++) {
    P->pos[i] =  P->L/((double)P->N-1)*i;
  }
  if(!(P->pres=malloc(P->N * sizeof(FLOAT)))){
    fprintf(stderr, "Problem with pressure allocation");
    exit(1);
  }
  init_to_zero(P->pres, P->N);
  if(!(P->rho=malloc(P->N * sizeof(FLOAT)))){
    fprintf(stderr, "Problem with pressure allocation");
    exit(1);
  }
  init_to_zero(P->rho, P->N);
  if(!(P->vel=malloc(P->N * sizeof(FLOAT)))){
    fprintf(stderr, "Problem with pressure allocation");
    exit(1);
  }
  init_to_zero(P->vel, P->N);
  if(!(P->ene=malloc(P->N * sizeof(FLOAT)))){
    fprintf(stderr, "Problem with pressure allocation");
    exit(1);
  }
  init_to_zero(P->ene, P->N);

  if(!(U->U_1=malloc(U->N * sizeof(FLOAT)))){
    fprintf(stderr, "Problem with U_1 allocation");
    exit(1);
  }
  init_to_zero(U->U_1, U->N);
  if(!(U->U_2=malloc(U->N * sizeof(FLOAT)))){
    fprintf(stderr, "Problem with U_1 allocation");
    exit(1);
  }
  init_to_zero(U->U_2, U->N);
  if(!(U->U_3=malloc(U->N * sizeof(FLOAT)))){
    fprintf(stderr, "Problem with U_1 allocation");
    exit(1);
  }
  init_to_zero(U->U_3, U->N);

  if(!(U->U_1_star=malloc(U->N * sizeof(FLOAT)))){
    fprintf(stderr, "Problem with U_1 allocation");
    exit(1);
  }
  init_to_zero(U->U_1_star, U->N);
  if(!(U->U_2_star=malloc(U->N * sizeof(FLOAT)))){
    fprintf(stderr, "Problem with U_1 allocation");
    exit(1);
  }
  init_to_zero(U->U_2_star, U->N);
  if(!(U->U_3_star=malloc(U->N * sizeof(FLOAT)))){
    fprintf(stderr, "Problem with U_1 allocation");
    exit(1);
  }
  init_to_zero(U->U_3_star, U->N);

  if(!(U->U_1_new=malloc(U->N * sizeof(FLOAT)))){
    fprintf(stderr, "Problem with U_1 allocation");
    exit(1);
  }
  init_to_zero(U->U_1_new, U->N);
  if(!(U->U_2_new=malloc(U->N * sizeof(FLOAT)))){
    fprintf(stderr, "Problem with U_1 allocation");
    exit(1);
  }
  init_to_zero(U->U_2_new, U->N);
  if(!(U->U_3_new=malloc(U->N * sizeof(FLOAT)))){
    fprintf(stderr, "Problem with U_1 allocation");
    exit(1);
  }
  init_to_zero(U->U_3_new, U->N);

  if(!(F->F_1=malloc(F->N * sizeof(FLOAT)))){
    fprintf(stderr, "Problem with F_1_X allocation");
    exit(1);
  }
  init_to_zero(F->F_1, F->N);
  if(!(F->F_2=malloc(F->N * sizeof(FLOAT)))){
    fprintf(stderr, "Problem with F_1_X allocation");
    exit(1);
  }
  init_to_zero(F->F_2, F->N);
  if(!(F->F_3=malloc(F->N * sizeof(FLOAT)))){
    fprintf(stderr, "Problem with F_1_X allocation");
    exit(1);
  }
  init_to_zero(F->F_3, F->N);
  if(!(F->F_1_star=malloc(F->N * sizeof(FLOAT)))){
    fprintf(stderr, "Problem with F_1_X allocation");
    exit(1);
  }
  init_to_zero(F->F_1_star, F->N);
  if(!(F->F_2_star=malloc(F->N * sizeof(FLOAT)))){
    fprintf(stderr, "Problem with F_1_X allocation");
    exit(1);
  }
  init_to_zero(F->F_2_star, F->N);
  if(!(F->F_3_star=malloc(F->N * sizeof(FLOAT)))){
    fprintf(stderr, "Problem with F_1_X allocation");
    exit(1);
  }
  init_to_zero(F->F_3_star, F->N);
}
void initial_conditions(physics_grid *P, U_grid *U, F_grid *F){
  int i;
  for (i = 0; i < P->N; i++) {
    P->vel[i] = V0;
    if(P->pos[i] <= X0){
      P->rho[i] = RL;
      P->pres[i] = PL;
    }
    else if(P->pos[i]> X0){
      P->rho[i] = RR;
      P->pres[i] = PR;
    }

    P->ene[i] = energy(P->pres[i],P->rho[i],P->vel[i]);

    U->U_1[i] = P->rho[i];
    U->U_2[i] = P->rho[i]*P->vel[i];
    U->U_3[i] = P->ene[i];
    F->F_1[i] = P->vel[i]*P->rho[i];
    F->F_2[i] = P->rho[i]*pow(P->vel[i],2.0) + P->pres[i];
    F->F_3[i] = P->vel[i]*(P->ene[i] + P->pres[i]);
  }
}
void fromU2F(F_grid *F, U_grid *U,int cases){
  int i;
  double v,p;
  // 1 == star
  if (cases == 1) {
    for (i = 0; i < F->N; i++) {
      v = U->U_2_star[i]/U->U_1_star[i];
      p = presion(U->U_3_star[i],U->U_1_star[i],v);
      F->F_1_star[i] = U->U_2_star[i];
      F->F_2_star[i] = U->U_1_star[i]*pow(v,2.0) + p;
      F->F_3_star[i] = v*(U->U_3_star[i] + p);
    }
  }
  // 0 == star
  else if (cases == 0){
    for (i = 0; i < F->N; i++) {
      v = U->U_2[i]/U->U_1[i];
      p = presion(U->U_3[i],U->U_1[i],v);
      F->F_1[i] = U->U_2[i];
      F->F_2[i] = U->U_1[i]*pow(v,2.0) + p;
      F->F_3[i] = v*(U->U_3[i] + p);
    }
  }
}
void U_star_update_mac(U_grid *U, F_grid *F,double alpha){
  U->U_1_star = U->U_1;
  U->U_2_star = U->U_2;
  U->U_3_star = U->U_3;
  int i;
  for ( i = 1; i < U->N-1; i++) {
    U->U_1_star[i]=U->U_1[i]-alpha*(F->F_1[i+1]-F->F_1[i]);
    U->U_2_star[i]=U->U_2[i]-alpha*(F->F_2[i+1]-F->F_2[i]);
    U->U_3_star[i]=U->U_3[i]-alpha*(F->F_3[i+1]-F->F_3[i]);
  }
}
void U_update_mac(U_grid *U, F_grid *F,double alpha){
  int i;
  U->U_1_new=U->U_1;
  U->U_2_new=U->U_2;
  U->U_3_new=U->U_3;
  for( i = 1; i<U->N;i++){
    U->U_1_new[i]=0.5*(U->U_1[i]+U->U_1_star[i]-alpha*(F->F_1_star[i]-F->F_1_star[i-1]));
    U->U_2_new[i]=0.5*(U->U_2[i]+U->U_2_star[i]-alpha*(F->F_2_star[i]-F->F_2_star[i-1]));
    U->U_3_new[i]=0.5*(U->U_3[i]+U->U_3_star[i]-alpha*(F->F_3_star[i]-F->F_3_star[i-1]));
  }
}
void U_star_update(U_grid *U, F_grid *F,double alpha){
  U->U_1_star = U->U_1;
  U->U_2_star = U->U_2;
  U->U_3_star = U->U_3;
  int i;
  for ( i = 1; i < U->N-1; i++) {
    U->U_1_star[i] = 0.5*((U->U_1[i+1]+ U->U_1[i]) - alpha*(F->F_1[i+1]-F->F_1[i]));
    U->U_2_star[i] = 0.5*((U->U_2[i+1]+ U->U_2[i]) - alpha*(F->F_2[i+1]-F->F_2[i]));
    U->U_3_star[i] = 0.5*((U->U_3[i+1]+ U->U_3[i]) - alpha*(F->F_3[i+1]-F->F_3[i]));
  }
}
void U_update(U_grid *U, F_grid *F,double alpha){
  int i;
  U->U_1_new=U->U_1;
  U->U_2_new=U->U_2;
  U->U_3_new=U->U_3;
  for( i = 1; i<U->N;i++){
    U->U_1_new[i] = (U->U_1[i] - alpha*(F->F_1_star[i]-F->F_1_star[i-1]));
    U->U_2_new[i] = (U->U_2[i] - alpha*(F->F_2_star[i]-F->F_2_star[i-1]));
    U->U_3_new[i] = (U->U_3[i] - alpha*(F->F_3_star[i]-F->F_3_star[i-1]));
  }
}
void fromU2var(physics_grid *P, U_grid *U){
  int i;
  for (i = 0; i < U->N; i++){
    P->rho[i]=U->U_1[i];
    P->vel[i]=U->U_2[i]/U->U_1[i];
    P->pres[i]=presion(U->U_3[i], U->U_1[i], P->vel[i]);
    P->ene[i]=U->U_3[i];
  }
}
double max_velocity(physics_grid * P){
  int i;
  double vel;
  double max_vel = sqrt((GAMMA*fabs(P->pres[0])/P->rho[0])) + fabs(P->vel[0]);
  for (i = 1; i < P->N; i++) {
    vel = sqrt((GAMMA*fabs(P->pres[i])/P->rho[i])) + fabs(P->vel[i]);
    if( vel > max_vel) {
      max_vel = vel;
    }
  }
  return max_vel;
}
