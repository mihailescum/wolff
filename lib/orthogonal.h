
#pragma once

#include <stdlib.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <cmath>

#include "state.h"
#include "types.h"

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
orthogonal_t <q, double> generate_rotation_uniform (gsl_rng *r, const state_t <orthogonal_t <q, double>, vector_t <q, double>> *s) {
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
orthogonal_t <q, double> generate_rotation_perturbation (gsl_rng *r, const state_t <orthogonal_t <q, double>, vector_t <q, double>> *s, double epsilon) {
  orthogonal_t <q, double> ptr;
  vector_t <q, double> tmp_v;
  ptr.is_reflection = true;

  tmp_v.x = (double *)malloc(q * sizeof(double));

  double tmp2 = 0;
  double M2 = 0;
  double tmpM = 0;

  for (q_t i = 0; i < q; i++) {
    tmp_v.x[i] = gsl_ran_ugaussian(r);
    M2 += pow(s->M.x[i], 2);
    tmpM += tmp_v.x[i] * s->M.x[i];
  }

  double v2 = 0;

  for (q_t i = 0; i < q; i++) {
    tmp_v.x[i] = (tmp_v.x[i] - tmpM * s->M.x[i] / M2) + epsilon * gsl_ran_ugaussian(r);
    v2 += tmp_v.x[i] * tmp_v.x[i];
  }

  double mag_v = sqrt(v2);

  for (q_t i = 0; i < q; i++) {
    tmp_v.x[i] /= mag_v;
  }

  vector_t <q, double> v = act(s->R, tmp_v);
  free(tmp_v.x);
  ptr.x = v.x;

  return ptr;
}

