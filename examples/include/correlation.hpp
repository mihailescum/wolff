
#pragma once

#include <wolff/types.h>
#include <wolff/state.hpp>

#include <fftw3.h>

template <class R_t, class X_t>
double correlation_length(const state_t <R_t, X_t>& s) {
  double total = 0;

#ifdef DIMENSION
  for (D_t j = 0; j < DIMENSION; j++) {
#else
  for (D_t j = 0; j < s.D; j++) {
#endif
    total += norm_squared(s.ReF[j]) + norm_squared(s.ImF[j]);
  }

  return total / s.D;
}

