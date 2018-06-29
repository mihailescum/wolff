
#pragma once

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "types.h"
#include "stack.h"

typedef struct {
  count_t n;
  double x;
  double x2;
  double x4;
  double m2;
  double m4;
  count_t W;
  double *xx;
  dll_t *x_window;
} meas_t;

typedef struct {
  uint64_t n;
  uint64_t W;
  double *OO;
  dll_t *Op;
  double O;
  double O2;
} autocorr_t;

void meas_update(meas_t *m, double x);

double meas_dx(const meas_t *m);

double meas_c(const meas_t *m);

double meas_dc(const meas_t *m);

void update_autocorr(autocorr_t *OO, double O);

double rho(const autocorr_t *o, uint64_t i);

void print_meas(const meas_t *m, const char *sym, FILE *outfile);
void print_vec_meas(q_t q, const meas_t **m, const char *sym, FILE *outfile);

void free_meas(meas_t *m);

meas_t *meas_initialize(count_t W);

double get_tau(const meas_t *m);

double Cxx(const meas_t *m, count_t t);

