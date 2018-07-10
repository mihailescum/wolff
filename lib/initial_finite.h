
#pragma once

#include <stdbool.h>

#include "types.h"
#include "dihedral.h"
#include "cluster_finite.h"

static char *finite_model_t_strings[] = {"ISING", "POTTS", "CLOCK", "DGM"};

typedef enum {
  ISING,
  POTTS,
  CLOCK,
  DGM
} finite_model_t;

state_finite_t *initial_finite_prepare_ising(D_t D, L_t L, double T, double *H);
state_finite_t *initial_finite_prepare_potts(D_t D, L_t L, q_t q, double T, double *H);
state_finite_t *initial_finite_prepare_clock(D_t D, L_t L, q_t q, double T, double *H);
state_finite_t *initial_finite_prepare_dgm(D_t D, L_t L, q_t q, double T, double *H);

void state_finite_free(state_finite_t *s);

double state_finite_energy(state_finite_t *s);

