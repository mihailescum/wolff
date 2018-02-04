
#pragma once

#include <assert.h>
#include <fftw3.h>
#include <float.h>
#include <getopt.h>
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

typedef struct {
  graph_t *g;
  q_t *spins;
  double T;
  double *J;
  double *H;
  double *J_probs;
  double *H_probs;
  double E;
  v_t *M;
  q_t q;
} ising_state_t;

v_t flip_cluster(ising_state_t *s, v_t v0, q_t s1, gsl_rng *r);

graph_t *graph_add_ext(const graph_t *g);

