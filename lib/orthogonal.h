
#pragma once

#include <stdlib.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <cmath>

#include "state.h"
#include "types.h"
#include "vector.h"

template <q_t q, class T>
struct orthogonal_t { bool is_reflection; T *x; };

template <q_t q, class T>
void init(orthogonal_t <q, T> *ptr) {
  ptr->is_reflection = false;
  ptr->x = (T *)calloc(q * q, sizeof(T));

  for (q_t i = 0; i < q; i++) {
    ptr->x[q * i + i] = (T)1;
  }
}

template <q_t q, class T>
orthogonal_t <q, T> copy (orthogonal_t <q, T> m) {
  orthogonal_t <q, T> m_copy;
  m_copy.is_reflection = m.is_reflection;

  q_t size;

  if (m.is_reflection) {
    size = q;
  } else {
    size = q * q;
  }

  m_copy.x = (T *)calloc(size, sizeof(T));

  for (q_t i = 0; i < size; i++) {
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

  if (m.is_reflection) {
    double prod = 0;
    for (q_t i = 0; i < q; i++) {
      prod += v.x[i] * m.x[i];
    }
    for (q_t i = 0; i < q; i++) {
      v_rot.x[i] = v.x[i] - 2 * prod * m.x[i];
    }
  } else {
    for (q_t i = 0; i < q; i++) {
      for (q_t j = 0; j < q; j++) {
        v_rot.x[i] += m.x[q * i + j] * v.x[j];
      }
    }
  }

  return v_rot;
}

template <q_t q, class T>
orthogonal_t <q, T> act (orthogonal_t <q, T> m1, orthogonal_t <q, T> m2) {
  orthogonal_t <q, T> m2_rot;

  m2_rot.is_reflection = false;
  m2_rot.x = (T *)calloc(q * q, sizeof(T));

  if (m1.is_reflection) {
    for (q_t i = 0; i < q; i++) {
      double akOki = 0;

      for (q_t k = 0; k < q; k++) {
        akOki += m1.x[k] * m2.x[q * k + i];
      }

      for (q_t j = 0; j < q; j++) {
        m2_rot.x[q * j + i] = m2.x[q * j + i] - 2 * akOki * m1.x[j];
      }
    }
  } else {
    for (q_t i = 0; i < q; i++) {
      for (q_t j = 0; j < q; j++) {
        for (q_t k = 0; k < q; k++) {
          m2_rot.x[i * q + j] += m1.x[i * q + j] * m2.x[j * q + k];
        }
      }
    }
  }

  return m2_rot;
}

template <q_t q, class T>
vector_t <q, T> act_inverse (orthogonal_t <q, T> m, vector_t <q, T> v) {
  if (m.is_reflection) {
    return act(m, v); // reflections are their own inverse
  } else {
    vector_t <q, T> v_rot;
    v_rot.x = (T *)calloc(q, sizeof(T));

    for (q_t i = 0; i < q; i++) {
      for (q_t j = 0; j < q; j++) {
        v_rot.x[i] += m.x[q * j + i] * v.x[j];
      }
    }

    return v_rot;
  }
}

template <q_t q, class T>
orthogonal_t <q, T> act_inverse (orthogonal_t <q, T> m1, orthogonal_t <q, T> m2) {
  if (m1.is_reflection) {
    return act(m1, m2); // reflections are their own inverse
  } else {
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
}

template <q_t q>
orthogonal_t <q, double> generate_rotation_uniform (gsl_rng *r, vector_t <q, double> v) {
  orthogonal_t <q, double> ptr;
  ptr.is_reflection = true;
  ptr.x = (double *)calloc(q, sizeof(double));

  double v2 = 0;

  for (q_t i = 0; i < q; i++) {
    ptr.x[i] = gsl_ran_ugaussian(r);
    v2 += ptr.x[i] * ptr.x[i];
  }

  double mag_v = sqrt(v2);

  for (q_t i = 0; i < q; i++) {
    ptr.x[i] /= mag_v;
  }

  return ptr;
}

template <q_t q>
orthogonal_t <q, double> generate_rotation_perturbation (gsl_rng *r, vector_t <q, double> v0, double epsilon, unsigned int n) {
  orthogonal_t <q, double> m;
  m.is_reflection = true;
  m.x = (double *)malloc(q * sizeof(double));

  vector_t <q, double> v;

  if (n > 1) {
    unsigned int rotation = gsl_rng_uniform_int(r, n);
    v.x = (double *)malloc(q * sizeof(double));

    double cosr = cos(2 * M_PI * rotation / (double)n / 2.0);
    double sinr = sin(2 * M_PI * rotation / (double)n / 2.0);

    v.x[0] = v0.x[0] * cosr - v0.x[1] * sinr;
    v.x[1] = v0.x[1] * cosr + v0.x[0] * sinr;

    for (q_t i = 2; i < q; i++) {
      v.x[i] = v0.x[i];
    }
  } else {
    v.x = v0.x;
  }

  double m2 = 0;
  double m_dot_v = 0;

  for (q_t i = 0; i < q; i++) {
    m.x[i] = gsl_ran_ugaussian(r);
    m_dot_v += m.x[i] * v.x[i];
  }

  double v2 = 0;
  double factor = epsilon * gsl_ran_ugaussian(r);

  for (q_t i = 0; i < q; i++) {
    m.x[i] = m.x[i] - m_dot_v * v.x[i] + factor * v.x[i];
    v2 += pow(m.x[i], 2);
  }

  double mag_v = sqrt(v2);

  for (q_t i = 0; i < q; i++) {
    m.x[i] /= mag_v;
  }

  if (n > 1) {
    free(v.x);
  }

  return m;
}

