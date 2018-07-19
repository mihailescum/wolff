
#pragma once

#include "types.h"
#include "state.h"

#include <fftw3.h>

template <class R_t, class X_t>
double correlation_length(const state_t <R_t, X_t> *s, fftw_plan plan, double *in, double *out) {
  for (v_t i = 0; i < s->nv; i++) {
    in[i] = correlation_component(s->spins[i]);
  }

  fftw_execute(plan);

  double length = pow(out[0], 2);

  return length;
}

