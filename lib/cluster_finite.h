
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
  graph_t *g;
  q_t q;
  R_t n_transformations;
  q_t *transformations;
  double T;
  double *J;
  double *H;
  double *J_probs;
  double *H_probs;
  q_t *spins;
  q_t *R;
  double E;
  v_t *M;
} state_finite_t;

v_t flip_cluster_finite(state_finite_t *s, v_t v0, q_t rot, gsl_rng *r);

