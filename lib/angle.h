#pragma once

#include "types.h"

#include <cmath>
#include "vector.h"

class angle_t {
  public:
    double x;

    typedef vector_t<2, double> M_t;
    typedef vector_t<2, double> F_t;

    angle_t() : x(0) {}
    angle_t(double x) : x(x) {}

    inline vector_t<2, double> operator*(v_t a) const {
      vector_t<2, double>M;
      M[0] = a * cos(x);
      M[1] = a * sin(x);

      return M;
    }

    inline vector_t<2, double> operator*(double a) const {
      vector_t<2, double>M;
      M[0] = a * cos(x);
      M[1] = a * sin(x);

      return M;
    }
};

inline vector_t<2,double>& operator+=(vector_t<2,double>& M, const angle_t& theta) {
  M[0] += cos(theta.x);
  M[1] += sin(theta.x);

  return M;
}

inline vector_t<2,double>& operator-=(vector_t<2,double>& M, const angle_t& theta) {
  M[0] -= cos(theta.x);
  M[1] -= sin(theta.x);

  return M;
}

