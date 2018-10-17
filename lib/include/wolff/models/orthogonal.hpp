
#pragma once

#include <random>
#include <cmath>

#include <wolff/types.h>
#include "vector.hpp"

template <q_t q, class T>
class orthogonal_t : public std::array<std::array<T, q>, q> {
  public :
  bool is_reflection;

  orthogonal_t() : is_reflection(false) {
    for (q_t i = 0; i < q; i++) {
      (*this)[i].fill(0);
      (*this)[i][i] = (T)1;
    }
  }

  vector_t<q, T> act(const vector_t <q, T>& v) const {
    vector_t <q, T> v_rot;
    v_rot.fill(0);

    if (is_reflection) {
      double prod = 0;
      for (q_t i = 0; i < q; i++) {
        prod += v[i] * (*this)[0][i];
      }
      for (q_t i = 0; i < q; i++) {
        v_rot[i] = v[i] - 2 * prod * (*this)[0][i];
      }
    } else {
      for (q_t i = 0; i < q; i++) {
        for (q_t j = 0; j < q; j++) {
          v_rot[i] += (*this)[i][j] * v[j];
        }
      }
    }

    return v_rot;
  }

  orthogonal_t<q, T> act(const orthogonal_t <q, T>& m) const {
    orthogonal_t <q, T> m_rot;

    m_rot.is_reflection = false;

    if (is_reflection) {
      for (q_t i = 0; i < q; i++) {
        double akOki = 0;

        for (q_t k = 0; k < q; k++) {
          akOki += (*this)[0][k] * m[k][i];
        }

        for (q_t j = 0; j < q; j++) {
          m_rot[j][i] = m[j][i] - 2 * akOki * (*this)[0][j];
        }
      }
    } else {
      for (q_t i = 0; i < q; i++) {
        m_rot[i].fill(0);
        for (q_t j = 0; j < q; j++) {
          for (q_t k = 0; k < q; k++) {
            m_rot[i][j] += (*this)[i][j] * m[j][k];
          }
        }
      }
    }

    return m_rot;
  }

  vector_t <q, T> act_inverse(const vector_t <q, T>& v) const {
    if (is_reflection) {
      return this->act(v); // reflections are their own inverse
    } else {
      vector_t <q, T> v_rot;
      v_rot.fill(0);

      for (q_t i = 0; i < q; i++) {
        for (q_t j = 0; j < q; j++) {
          v_rot[i] += (*this)[j][i] * v[j];
        }
      }

      return v_rot;
    }
  }

  vector_t <q, T> act_inverse(const orthogonal_t <q, T>& m) const {
    if (is_reflection) {
      return this->act(m); // reflections are their own inverse
    } else {
      orthogonal_t <q, T> m_rot;
      m_rot.is_reflection = false;

      for (q_t i = 0; i < q; i++) {
        m_rot[i].fill(0);
        for (q_t j = 0; j < q; j++) {
          for (q_t k = 0; k < q; k++) {
            m_rot[i][j] += (*this)[j][i] * m[j][k];
          }
        }
      }

      return m_rot;
    }
  }

};

template <q_t q>
orthogonal_t <q, double> generate_rotation_uniform (std::mt19937& r, const vector_t <q, double>& v) {
  std::normal_distribution<double> dist(0.0,1.0);
  orthogonal_t <q, double> ptr;
  ptr.is_reflection = true;

  double v2 = 0;

  for (q_t i = 0; i < q; i++) {
    ptr[0][i] = dist(r);
    v2 += ptr[0][i] * ptr[0][i];
  }

  double mag_v = sqrt(v2);

  for (q_t i = 0; i < q; i++) {
    ptr[0][i] /= mag_v;
  }

  return ptr;
}

template <q_t q>
orthogonal_t <q, double> generate_rotation_perturbation (std::mt19937& r, const vector_t <q, double>& v0, double epsilon, unsigned int n) {
  std::normal_distribution<double> dist(0.0,1.0);
  orthogonal_t <q, double> m;
  m.is_reflection = true;

  vector_t <q, double> v;

  if (n > 1) {
    std::uniform_int_distribution<unsigned int> udist(0, n);
    unsigned int rotation = udist(r);

    double cosr = cos(2 * M_PI * rotation / (double)n / 2.0);
    double sinr = sin(2 * M_PI * rotation / (double)n / 2.0);

    v[0] = v0[0] * cosr - v0[1] * sinr;
    v[1] = v0[1] * cosr + v0[0] * sinr;

    for (q_t i = 2; i < q; i++) {
      v[i] = v0[i];
    }
  } else {
    v = v0;
  }

  double m_dot_v = 0;

  for (q_t i = 0; i < q; i++) {
    m[0][i] = dist(r); // create a random vector
    m_dot_v += m[0][i] * v[i];
  }

  double v2 = 0;

  for (q_t i = 0; i < q; i++) {
    m[0][i] = m[0][i] - m_dot_v * v[i]; // find the component orthogonal to v
    v2 += pow(m[0][i], 2);
  }

  double mag_v = sqrt(v2);

  for (q_t i = 0; i < q; i++) {
    m[0][i] /= mag_v; // normalize
  }

  v2 = 0;

  double factor = epsilon * dist(r);

  for (q_t i = 0; i < q; i++) {
    m[0][i] += factor * v[i]; // perturb orthogonal vector in original direction
    v2 += pow(m[0][i], 2);
  }

  mag_v = sqrt(v2);

  for (q_t i = 0; i < q; i++) {
    m[0][i] /= mag_v; // normalize
  }

  return m;
}

