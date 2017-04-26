#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "struct.h"
#include "init.h"
#include "io.h"

int main(int argc, char const *argv[]) {
  double t=0.0;
  int i;
  double delta_t=0.0001;
  double T = 0.1;
  double alpha;
  physics_grid * P_state;
  U_grid * U_state;
  F_grid * F_state;

  P_state = create_physics_grid();
  U_state = create_U_grid();
  F_state = create_F_grid();

  init_problem(P_state, U_state, F_state);
  initial_conditions(P_state, U_state, F_state);

  while (t<T) {
    delta_t = (0.5*DELTA_X)/max_velocity(P_state);
    alpha = delta_t/P_state->delta_x;
    fromU2F(F_state,U_state,NO_STAR);
    U_star_update(U_state,F_state,alpha);
    fromU2F(F_state,U_state,STAR);
    U_update(U_state,F_state,alpha);
    U_state->U_1 = U_state->U_1_new;
    U_state->U_2 = U_state->U_2_new;
    U_state->U_3 = U_state->U_3_new;
    t = t + delta_t;
  }
  fromU2var(P_state,U_state);
  /* ------------------------------- PRINT ---------------------------------- */
  print_physics_grid(P_state);
  return 0;
}
