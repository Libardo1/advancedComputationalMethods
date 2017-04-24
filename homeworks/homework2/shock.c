#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double gam =1.4;
double deltat = 0.001;
double deltax = 0.01;
double L = 1.0;
double T = 0.01;
int N;
double alpha;

/* -----------------------------  Declarations -------------------------------*/
void init_to_zero(double *p, int n_points);
void init_physc(double * pos,double *rho,double * pres,double * vel);
void init_grids(double *H1, double *H2, double *H3);
double energy(double p,double r,double v);
double presion(double e, double r,double v);
void init(double *r, double *p, double *v,double*e, double*x, double *u1, double *u2, double *u3, double *f1,double *f2, double *uf3);
void fromU2F(double*f1, double*f2, double*f3,double*u1,double*u2,double*u3);
/* Lax Wendorf functions*/
void Ustar_update(double *u1_star,double *u2_star,double *u3_star,double *u1,double *u2,double *u3,double *f1,double *f2,double *f3,double r);
void U_update(double *u1_new,double *u2_new,double *u3_new,double *u1,double *u2,double *u3, double *u1_star,double *u2_star,double *u3_star,double *f1_star,double *f2_star,double *f3_star,double r);
void fromU2var(double*r,double*v,double*p,double*e,double*u1,double*u2,double*u3);

/* -----------------------------      Main     -------------------------------*/

int main(int argc, char const *argv[]) {
  int i;
  double t;

  N = (int) L/deltax;

  double * U1;
  double * U2;
  double * U3;
  double * U1_star;
  double * U2_star;
  double * U3_star;
  double * U1_new;
  double * U2_new;
  double * U3_new;
  double * F1;
  double * F2;
  double * F3;
  double * F1_star;
  double * F2_star;
  double * F3_star;

  double * pos;
  double * rho;
  double * pres;
  double * vel;
  double * e;

  /* -----------------------------  malloc -------------------------------*/

  U1=malloc(N * sizeof(double));
  U2=malloc(N * sizeof(double));
  U3=malloc(N * sizeof(double));
  U1_star=malloc(N * sizeof(double));
  U2_star=malloc(N * sizeof(double));
  U3_star=malloc(N * sizeof(double));
  U1_new=malloc(N * sizeof(double));
  U2_new=malloc(N * sizeof(double));
  U3_new=malloc(N * sizeof(double));
  F1=malloc(N * sizeof(double));
  F2=malloc(N * sizeof(double));
  F3=malloc(N * sizeof(double));
  F1_star=malloc(N * sizeof(double));
  F2_star=malloc(N * sizeof(double));
  F3_star=malloc(N * sizeof(double));
  pos=malloc(N * sizeof(double));
  rho=malloc(N * sizeof(double));
  pres=malloc(N * sizeof(double));
  vel=malloc(N * sizeof(double));
  e=malloc(N * sizeof(double));

  init_to_zero(U1,N);
  init_to_zero(U2,N);
  init_to_zero(U3,N);
  init_to_zero(U1_star,N);
  init_to_zero(U2_star,N);
  init_to_zero(U3_star,N);
  init_to_zero(U1_new,N);
  init_to_zero(U2_new,N);
  init_to_zero(U3_new,N);
  init_to_zero(F1,N);
  init_to_zero(F2,N);
  init_to_zero(F3,N);
  init_to_zero(F1_star,N);
  init_to_zero(F2_star,N);
  init_to_zero(F3_star,N);
  init_to_zero(pos,N);
  init_to_zero(rho,N);
  init_to_zero(pres,N);
  init_to_zero(e,N);
  for ( i = 0; i < N; i++) {
    pos[i] =  L/((double)N-1)*i;
  }

  init(rho,pres,vel,e,pos,U1,U2,U3,F1,F2,F3);
  for (i = 0; i < N; i++) {
    //printf("%f %f %f %f %f \n", pos[i],rho[i],vel[i],pres[i],e[i]);
  }
  while (t<deltat) {
    alpha = deltat/deltax;
    Ustar_update(U1_star,U2_star,U3_star,U1,U2,U3,F1,F2,F3,alpha);
    for (i = 0; i < N; i++) {
      printf("%f %f %f \n", U1_star[i],U2_star[i],U3_star[i]);
    }

    fromU2F(F1_star,F2_star,F3_star,U1_star,U2_star,U3_star);
    U_update(U1_new,U2_new,U3_new,U1,U2,U3,U1_star,U2_star,U3_star,F1_star,F2_star,F3_star,alpha);
    U1 = U1_new;
    U2 = U2_new;
    U3 = U3_new;
    t = t + deltat;
  }
  fromU2var(rho,vel,pres,e,U1,U2,U3);
  for (i = 0; i < N; i++) {
    printf("%f %f %f %f %f \n", pos[i],rho[i],vel[i],pres[i],e[i]);
  }
  return 0;
}

void init_to_zero(double * p, int n_points){
  int i;
  for(i=0;i<n_points;i++){
    p[i] = 0.0;
  }
}
void init_physc(double * pos,double *rho,double * pres,double * vel){
  int i;
  if(!(pos=malloc(N * sizeof(double)))){
     fprintf(stderr, "Problem with pressure allocation");
     exit(1);
  }

  if(!(rho=malloc(N * sizeof(double)))){
   fprintf(stderr, "Problem with pressure allocation");
   exit(1);
  }

  if(!(pres=malloc(N * sizeof(double)))){
     fprintf(stderr, "Problem with pressure allocation");
     exit(1);
  }

  if(!(vel=malloc(N * sizeof(double)))){
   fprintf(stderr, "Problem with pressure allocation");
   exit(1);
  }

  for ( i = 0; i < N; i++) {
    pos[i] =  L/((double)N-1)*i;
  }
  init_to_zero(rho,N);
  init_to_zero(pres,N);
  init_to_zero(vel,N);

}
void init_grids(double* H1,double* H2,double* H3){
  if(!(H1=malloc(N * sizeof(double)))){
     fprintf(stderr, "Problem with pressure allocation");
     exit(1);
  }
  if(!(H2=malloc(N * sizeof(double)))){
     fprintf(stderr, "Problem with pressure allocation");
     exit(1);
  }
  if(!(H3=malloc(N * sizeof(double)))){
     fprintf(stderr, "Problem with pressure allocation");
     exit(1);
  }
  init_to_zero(H1,N);
  init_to_zero(H2,N);
  init_to_zero(H3,N);
}
void init(double *r, double *p, double *v,double*e, double*x, double *u1, double *u2, double *u3, double *f1,double *f2, double *f3){
  int i;

  for (i = 0; i < N; i++) {
    v[i] = 0.0;

    if(x[i]<= 0.5){
      r[i] = 1.0;
      p[i] = 1.0;
    }
    else if(x[i]>0.5){
      r[i] = 0.125;
      p[i] = 0.10;
    }

    e[i] = energy(p[i],r[i],v[i]);
    u1[i] = r[i];
    u2[i] = r[i]*v[i];
    u3[i] = e[i];
    f1[i] = v[i]*r[i];
    f2[i] = r[i]*pow(v[i],2.0) + p[i];
    f3[i] = v[i]*e[i] + p[i];
  }

}

double energy(double p,double r,double v){
  return (p/(gam-1)) + (r*pow(v,2.0))/2.0;
}
double presion(double e, double r,double v){
    return (gam-1)*(e - r *pow(v,2.0))/2.0;
}
void fromU2F(double*f1, double*f2, double*f3,double*u1,double*u2,double*u3){
  int i;
  double v, p;

  for (i = 0; i < N; i++) {
    v = u2[i]/u1[i];
    p = presion(u3[i],u1[i],v);
    f1[i] = u2[i];
    f2[i] = u1[i]*pow(v,2.0) + p;
    f3[i] = v*u3[i] + p;
  }

}
/* Lax Wendorf functions*/
void Ustar_update(double *u1_star,double *u2_star,double *u3_star,double *u1,double *u2,double *u3,double *f1,double *f2,double *f3,double alpha) {
  u1_star = u1;
  u2_star = u2;
  u3_star = u3;
  int i;

  for ( i = 1; i < N-1; i++) {
    u1_star[i]=u1[i]-alpha*(f1[i+1]-f1[i]);
    u2_star[i]=u2[i]-alpha*(f2[i+1]-f2[i]);
    u3_star[i]=u3[i]-alpha*(f3[i+1]-f3[i]);
  }
}
void U_update(double *u1_new,double *u2_new,double *u3_new,double *u1,double *u2,double *u3, double *u1_star,double *u2_star,double *u3_star,double *f1_star,double *f2_star,double *f3_star,double alpha){
  u1_new=u1;
  u2_new=u2;
  u3_new=u3;
  int i;
  for( i = 1; i<N;i++){
    u1_new[i]=0.5*(u1[i]+u1_star[i]-alpha*(f1_star[i]-f1_star[i-1]));
    u2_new[i]=0.5*(u2[i]+u2_star[i]-alpha*(f2_star[i]-f2_star[i-1]));
    u3_new[i]=0.5*(u3[i]+u3_star[i]-alpha*(f3_star[i]-f3_star[i-1]));
  }

}
void fromU2var(double*r,double*v,double*p,double*e,double*u1,double*u2,double*u3){
  int i;
  for ( i = 0; i < N; i++) {
    r[i]=u1[i];
    v[i]=u2[i]/u1[i];
    p[i]=presion(u3[i], u1[i], v[i]);
    e[i]=u3[i];
  }
}
