
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

#include <jst/graph.h>
#include <jst/rand.h>

#include "queue.h"

typedef struct {
  graph_t *g;
  bool *spins;
  int32_t M;
  double H;
} ising_state_t;

typedef struct {
  uint32_t nv;
  int32_t dJb;
  int32_t dHb;
} cluster_t;

int32_t sign(double x);

cluster_t *flip_cluster(const graph_t *g, const double *ps, bool *x,
                        gsl_rng *r);

graph_t *graph_add_ext(const graph_t *g);

uint32_t wolff_step(double T, double H, ising_state_t *s, gsl_rng *r,
                    double *ps);
