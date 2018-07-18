
#pragma once

#include "types.h"
#include "state.h"

#include <fftw3.h>

template <class R_t, class X_t>
double correlation_length(const state_t <R_t, X_t> *s) {
  double *data = (double *)fftw_malloc(s->nv * sizeof(double));
  int rank = s->D;
  int *n = (int *)malloc(rank * sizeof(int));
  fftw_r2r_kind *kind = (fftw_r2r_kind *)malloc(rank * sizeof(fftw_r2r_kind));
  for (D_t i = 0; i < rank; i++) {
    n[i] = s->L;
    kind[i] = FFTW_R2HC;
  }
  fftw_plan plan = fftw_plan_r2r(rank, n, data, data, kind, 0);

  for (v_t i = 0; i < s->nv; i++) {
    data[i] = correlation_component(s->spins[i]);
  }

  fftw_execute(plan);

  double length = pow(data[0], 2);

  fftw_destroy_plan(plan);
  fftw_free(data);
  free(n);
  free(kind);

  return length;
}

