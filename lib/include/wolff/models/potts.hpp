
#pragma once

#include <cmath>

#include "../types.h"
#include "vector.hpp"

namespace wolff {

template <q_t q>
class potts_t {
  public:
    q_t x;

    potts_t() : x(0) {}
    potts_t(q_t x) : x(x) {}

    typedef vector_t<q, int> M_t;
    typedef vector_t<q, double> F_t;

    inline vector_t<q, int> operator*(v_t a) const {
      vector_t<q, int> result;
      result.fill(0);
      result[x] = (int)a;

      return result;
    }

    inline vector_t<q, double> operator*(double a) const {
      vector_t<q, double> result;
      result.fill(0.0);
      result[x] = a;

      return result;
    }

    inline vector_t<q, int> operator-(const potts_t<q> &s) const {
      vector_t<q, int> result;
      result.fill(0);

      result[x]++;
      result[s.x]--;

      return result;
    }

    q_t enumerate() const {
      return x;
    }
};

template<q_t q>
inline typename potts_t<q>::M_t operator*(v_t a, const potts_t<q>& s) {
  return s * a;
}

template<q_t q>
inline typename potts_t<q>::F_t operator*(double a, const potts_t<q>& s) {
  return s * a;
}

}

