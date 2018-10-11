
#pragma once

#include <stdlib.h>
#include <cmath>
#include <array>

#include <wolff/types.h>

template <q_t q, class T>
class vector_t : public std::array<T, q> {
  public: 

    // M_t needs to hold the sum of nv spins
    typedef vector_t <q, T> M_t;

    // F_t needs to hold the double-weighted sum of spins
    typedef vector_t <q, double> F_t;

    vector_t() {
      this->fill((T)0);
      (*this)[1] = (T)1;
    }

    vector_t(const T *x) {
      for (q_t i = 0; i < q; i++) {
        (*this)[i] = x[i];
      }
    }

    template <class U>
    inline vector_t<q, T>& operator+=(const vector_t<q, U> &v) {
      for (q_t i = 0; i < q; i++) {
        (*this)[i] += (U)v[i];
      }
      return *this;
    }

    template <class U>
    inline vector_t<q, T>& operator-=(const vector_t<q, U> &v) {
      for (q_t i = 0; i < q; i++) {
        (*this)[i] -= (U)v[i];
      }
      return *this;
    }

    inline vector_t<q, T> operator*(v_t x) const {
      vector_t<q, T> result;
      for (q_t i = 0; i < q; i++) {
        result[i] = x * (*this)[i];
      }

      return result;
    }

    inline vector_t<q, double> operator*(double x) const {
      vector_t<q, double> result;
      for (q_t i = 0; i < q; i++) {
        result[i] = x * (*this)[i];
      }

      return result;
    }

    inline vector_t<q, T> operator-(const vector_t<q, T>& v) const {
      vector_t<q, T> diff = *this;
      diff -= v;
      return diff;
    }
};


template<q_t q, class T>
double norm_squared(vector_t<q, T> v) {
  double tmp = 0;
  for (T &x : v) {
    tmp += pow(x, 2);
  }

  return tmp;
}

template <q_t q, class T>
void write_magnetization(vector_t <q, T> M, FILE *outfile) {
  for (q_t i = 0; i < q; i++) {
    fwrite(&(M[i]), sizeof(T), q, outfile);
  }
}

// below functions and definitions are unnecessary for wolff.h but useful.

template <q_t q> // save some space and don't write whole doubles
void write_magnetization(vector_t <q, double> M, FILE *outfile) {
  for (q_t i = 0; i < q; i++) {
    float M_tmp = (float)M[i];
    fwrite(&M_tmp, sizeof(float), 1, outfile);
  }
}

template <q_t q, class T>
T dot(const vector_t <q, T>& v1, const vector_t <q, T>& v2) {
  T prod = 0;

  for (q_t i = 0; i < q; i++) {
    prod += v1[i] * v2[i];
  }

  return prod;
}

template <q_t q, class T>
double H_vector(const vector_t <q, T>& v1, T *H) {
  vector_t <q, T> H_vec(H);
  return (double)(dot <q, T> (v1, H_vec));
}

char const *ON_strings[] = {"TRIVIAL", "ISING", "PLANAR", "HEISENBERG"};

