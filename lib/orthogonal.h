
#pragma once

#include <stdlib.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <cmath>

#include "types.h"

template <q_t q, class T>
struct orthogonal_t { T *x; };

template <q_t q, class T>
void init(orthogonal_t <q, T> *ptr) {
  ptr->x = (T *)calloc(q * q, sizeof(T));

  for (q_t i = 0; i < q; i++) {
    ptr->x[q * i + i] = (T)1;
  }
}

template <q_t q, class T>
orthogonal_t <q, T> copy (orthogonal_t <q, T> m) {
  orthogonal_t <q, T> m_copy;
  m_copy.x = (T *)calloc(q * q, sizeof(T));

  for (q_t i = 0; i < q * q; i++) {
    m_copy.x[i] = m.x[i];
  }

  return m_copy;
}

template <q_t q, class T>
void free_spin (orthogonal_t <q, T> m) {
  free(m.x);
}

template <q_t q, class T>
vector_t <q, T> act (orthogonal_t <q, T> m, vector_t <q, T> v) {
  vector_t <q, T> v_rot;
  v_rot.x = (T *)calloc(q, sizeof(T));

  for (q_t i = 0; i < q; i++) {
    for (q_t j = 0; j < q; j++) {
      v_rot.x[i] += m.x[q * i + j] * v.x[j];
    }
  }

  return v_rot;
}

template <q_t q, class T>
orthogonal_t <q, T> act (orthogonal_t <q, T> m1, orthogonal_t <q, T> m2) {
  orthogonal_t <q, T> m2_rot;
  m2_rot.x = (T *)calloc(q * q, sizeof(T));

  for (q_t i = 0; i < q; i++) {
    for (q_t j = 0; j < q; j++) {
      for (q_t k = 0; k < q; k++) {
        m2_rot.x[i * q + j] += m1.x[i * q + j] * m2.x[j * q + k];
      }
    }
  }

  return m2_rot;
}

template <q_t q, class T>
vector_t <q, T> act_inverse (orthogonal_t <q, T> m, vector_t <q, T> v) {
  vector_t <q, T> v_rot;
  v_rot.x = (T *)calloc(q, sizeof(T));

  for (q_t i = 0; i < q; i++) {
    for (q_t j = 0; j < q; j++) {
      v_rot.x[i] += m.x[q * j + i] * v.x[j];
    }
  }

  return v_rot;
}

template <q_t q, class T>
orthogonal_t <q, T> act_inverse (orthogonal_t <q, T> m1, orthogonal_t <q, T> m2) {
  orthogonal_t <q, T> m2_rot;
  m2_rot.x = (T *)calloc(q * q, sizeof(T));

  for (q_t i = 0; i < q; i++) {
    for (q_t j = 0; j < q; j++) {
      for (q_t k = 0; k < q; k++) {
        m2_rot.x[i * q + j] += m1.x[j * q + i] * m2.x[j * q + k];
      }
    }
  }

  return m2_rot;
}

template <q_t q>
void generate_rotation (gsl_rng *r, orthogonal_t <q, double> *ptr) {
  double *v = (double *)malloc(q * sizeof(double));
  double v2 = 0;

  for (q_t i = 0; i < q; i++) {
    v[i] = gsl_ran_ugaussian(r);
    v2 += v[i] * v[i];
  }

  ptr->x = (double *)calloc(q * q, sizeof(double));
  
  for (q_t i = 0; i < q; i++) {
    ptr->x[q * i + i] = 1.0;
    for (q_t j = 0; j < q; j++) {
      ptr->x[q * i + j] -= 2 * v[i] * v[j] / v2;
    }
  }

  free(v);
}


