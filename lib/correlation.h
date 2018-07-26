
#pragma once

#include "types.h"
#include "state.h"

#include <fftw3.h>

template <class R_t, class X_t>
double correlation_length(const state_t <R_t, X_t>& s) {
  double total = 0;

  for (D_t j = 0; j < s.D; j++) {
    total += norm_squared(s.ReF[j]) + norm_squared(s.ImF[j]);
  }

  return total / s.D;
}

