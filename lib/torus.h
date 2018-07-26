
#pragma once

#include <cmath>
#include <array>
#include "types.h"

template <q_t n>
class torus_t : public std::array<double, n> {
  public:
    typedef std::array<double, n> M_t;
    typedef std::array<double, n> F_t;

    torus_t() {
      this->fill(0);
    }

    inline torus_t<n> operator*(v_t a) const {
      torus_t<n> x;
      for (q_t i = 0; i < n; i++) {
        x[i] = a * (*this)[i];
      }

      return x;
    }

    inline torus_t<n> operator*(double a) const {
      torus_t<n> x;
      for (q_t i = 0; i < n; i++) {
        x[i] = a * (*this)[i];
      }

      return x;
    }

    inline torus_t<n>& operator+=(const torus_t<n>& x) {
      for (q_t i = 0; i < n; i++) {
        (*this)[i] += x[i];
      }
    }

    inline torus_t<n>& operator-=(const torus_t<n>& x) {
      for (q_t i = 0; i < n; i++) {
        (*this)[i] -= x[i];
      }
    }
};

template <q_t n>
double norm_squared(const torus_t<n>& x) {
  double tmp = 0;
  for (const double& xi : x) {
    tmp += pow(xi, 2);
  }
  return tmp;
}

void write_magnetization(const torus_t<n>& x, FILE *outfile) {
  for (const double& xi : x) {
    float tmp_xi = (float)xi;
    fwrite(&tmp_xi, sizeof(float), 1, outfile);
  }
}

