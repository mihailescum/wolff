
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

typedef struct {
  graph_t *g;
  q_t *spins;
  double T;
  double *H;
  double T_prob;
  double *H_probs;
  double E;
  v_t *M;
  q_t q;
} ising_state_t;

typedef struct {
  uint64_t n;
  double x;
  double dx;
  double x2;
  double m2;
  double m4;
  double c;
  double dc;
} meas_t;

typedef struct {
  uint64_t n;
  uint64_t W;
  double *OO;
  dll_t *Op;
  double O;
  double O2;
} autocorr_t;

v_t flip_cluster(ising_state_t *s, v_t v0, q_t s1, gsl_rng *r);

graph_t *graph_add_ext(const graph_t *g);

void update_meas(meas_t *m, double x);

void update_autocorr(autocorr_t *OO, double O);

double rho(autocorr_t *o, uint64_t i);

