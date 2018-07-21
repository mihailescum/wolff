
#pragma once

#include "measurement.h"

#define POSSIBLE_MEASUREMENTS 5
const unsigned char measurement_energy        = 1 << 0;
const unsigned char measurement_clusterSize   = 1 << 1;
const unsigned char measurement_magnetization = 1 << 2;
const unsigned char measurement_fourierZero    = 1 << 3;

#ifdef __cplusplus

#include "state.h"
#include "correlation.h"
#include <functional>

template <class R_t, class X_t>
std::function <void(const state_t <R_t, X_t> *)> measurement_energy_file(FILE *file) {
  return [=](const state_t <R_t, X_t> *s) {
    float smaller_E = (float)s->E;
    fwrite(&smaller_E, sizeof(float), 1, file);
  };
}

template <class R_t, class X_t>
std::function <void(const state_t <R_t, X_t> *)> measurement_cluster_file(FILE *file) {
  return [=](const state_t <R_t, X_t> *s) {
    fwrite(&(s->last_cluster_size), sizeof(uint32_t), 1, file);
  };
}

template <class R_t, class X_t>
std::function <void(const state_t <R_t, X_t> *)> measurement_magnetization_file(FILE *file) {
  return [=](const state_t <R_t, X_t> *s) {
    write_magnetization(s->M, file);
  };
}

template <class R_t, class X_t>
std::function <void(const state_t <R_t, X_t> *)> measurement_fourier_file(FILE *file) {
  return [=](const state_t <R_t, X_t> *s) {
    float smaller_X = (float)correlation_length(s);
    fwrite(&smaller_X, sizeof(float), 1, file);
  };
}

template <class R_t, class X_t>
std::function <void(const state_t <R_t, X_t> *)> measurement_average_cluster(meas_t *x) {
  return [=](const state_t <R_t, X_t> *s) {
    meas_update(x, s->last_cluster_size);
  };
}

#endif


