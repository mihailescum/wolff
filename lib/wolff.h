
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

typedef enum {
  WOLFF,
  WOLFF_GHOST,
  METROPOLIS
} sim_t;

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
  bool hit_ghost;
  ll_t *spins;
} cluster_t;

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
  double *Op;
  double O;
  double O2;
} autocorr_t;

int8_t sign(double x);

cluster_t *flip_cluster(const graph_t *g, const double *ps, bool *x, bool stop_on_ghost,
                        gsl_rng *r);

graph_t *graph_add_ext(const graph_t *g);

uint32_t wolff_step(double T, double H, ising_state_t *s, sim_t sim, gsl_rng *r,
                    double *ps);

void update_meas(meas_t *m, double x);

void update_autocorr(autocorr_t *OO, double O);

double add_to_avg(double mx, double x, uint64_t n);

double rho(autocorr_t *o, uint64_t i);

