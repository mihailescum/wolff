
#pragma once

#include <assert.h>
#include <fftw3.h>
#include <float.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <inttypes.h>
#include <math.h>
#include <stdbool.h>
#include <string.h>
#include <sys/types.h>

#include "types.h"
#include "rand.h"
#include "stack.h"
#include "convex.h"
#include "graph.h"
#include "tree.h"
#include "measurement.h"
#include "symmetric.h"
#include "yule_walker.h"

typedef struct {
  D_t D;
  L_t L;
  v_t nv;
  v_t ne;
  graph_t *g;
  q_t q;
  R_t n_transformations;
  q_t *transformations;
  R_t n_involutions;
  R_t *involutions;
  R_t *transform_site_to_zero;
  q_t n_bond_types;
  q_t *bond_with_zero_type;
  double T;
  double *J;
  double *H;
  double *J_probs;
  double *H_probs;
  q_t *spins;
  q_t *R;
  v_t *B;
  v_t *M;
  double *ReF;
  double *ImF;
  double *precomputed_cos;
  double *precomputed_sin;
} state_finite_t;

v_t flip_cluster_finite(state_finite_t *s, v_t v0, R_t rot, gsl_rng *r);

