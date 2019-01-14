
#ifndef WOLFF_MODELS_POTTS_H
#define WOLFF_MODELS_POTTS_H

#include <cmath>

#include "vector.hpp"

namespace wolff {

  template <unsigned q>
  class potts_t {
    public:
      unsigned x;

      potts_t() : x(0) {}
      potts_t(unsigned x) : x(x) {}

      typedef vector_t<q, int> M_t;
      typedef vector_t<q, double> F_t;

      inline vector_t<q, int> operator*(unsigned a) const {
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

      unsigned enumerate() const {
        return x;
      }
  };

  template<unsigned q>
  inline typename potts_t<q>::M_t operator*(unsigned a, const potts_t<q>& s) {
    return s * a;
  }

  template<unsigned q>
  inline typename potts_t<q>::F_t operator*(double a, const potts_t<q>& s) {
    return s * a;
  }

}

#endif

