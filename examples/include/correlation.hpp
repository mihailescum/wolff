
#pragma once

#include <wolff/types.h>
#include <wolff/state.hpp>

#include <fftw3.h>

template <class X_t>
double correlation_length(const std::vector<typename X_t::F_t>& ReF, const std::vector<typename X_t::F_t>& ImF, D_t D) {
  double total = 0;

#ifdef DIMENSION
  for (D_t j = 0; j < DIMENSION; j++) {
#else
  for (D_t j = 0; j < D; j++) {
#endif
    total += norm_squared(ReF[j]) + norm_squared(ImF[j]);
  }

  return total / D;
}

