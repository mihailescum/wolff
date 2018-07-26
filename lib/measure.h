
#pragma once

#include "state.h"
#include "correlation.h"
#include <functional>

#define POSSIBLE_MEASUREMENTS 4
const unsigned char measurement_energy        = 1 << 0;
const unsigned char measurement_clusterSize   = 1 << 1;
const unsigned char measurement_magnetization = 1 << 2;
const unsigned char measurement_fourierZero    = 1 << 3;

char const *measurement_labels[] = {"E", "S", "M", "F"};

FILE **measure_setup_files(unsigned char flags, unsigned long timestamp) {
  FILE **files = (FILE **)calloc(POSSIBLE_MEASUREMENTS, sizeof(FILE *));

  for (uint8_t i = 0; i < POSSIBLE_MEASUREMENTS; i++) {
    if (flags & (1 << i)) {
      char *filename = (char *)malloc(255 * sizeof(char));
      sprintf(filename, "wolff_%lu_%s.dat", timestamp, measurement_labels[i]);
      files[i] = fopen(filename, "wb");
      free(filename);
    }
  }

  return files;
}

template <class R_t, class X_t>
std::function <void(const state_t <R_t, X_t>&)> measure_function_write_files(unsigned char flags, FILE **files, std::function <void(const state_t <R_t, X_t>&)> other_f) {
  return [=] (const state_t <R_t, X_t>& s) {
    if (flags & measurement_energy) {
      float smaller_E = (float)s.E;
      fwrite(&smaller_E, sizeof(float), 1, files[0]);
    }
    if (flags & measurement_clusterSize) {
      fwrite(&(s.last_cluster_size), sizeof(uint32_t), 1, files[1]);
    }
    if (flags & measurement_magnetization) {
      write_magnetization(s.M, files[2]);
    }
    if (flags & measurement_fourierZero) {
      float smaller_X = (float)correlation_length(s);
      fwrite(&smaller_X, sizeof(float), 1, files[3]);
    }

    other_f(s);
  };
}

void measure_free_files(unsigned char flags, FILE **files) {
  for (uint8_t i = 0; i < POSSIBLE_MEASUREMENTS; i++) {
    if (flags & (1 << i)) {
      fclose(files[i]);
    }
  }

  free(files);
}


