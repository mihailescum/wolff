
#pragma once

#include <cmath>
#include <array>
#include <iostream>

#include <wolff/types.h>

template <q_t q, class T>
class vector_t : public std::array<T, q> {
  public: 

    vector_t() {
      this->fill((T)0);
      (*this)[0] = (T)1;
    }

    vector_t(const T *x) {
      for (q_t i = 0; i < q; i++) {
        (*this)[i] = x[i];
      }
    }

    typedef vector_t <q, T> M_t;
    typedef vector_t <q, double> F_t;

    template <class U>
    inline vector_t<q, T>& operator+=(const vector_t<q, U> &v) {
      for (q_t i = 0; i < q; i++) {
        (*this)[i] += (T)v[i];
      }
      return *this;
    }

    template <class U>
    inline vector_t<q, T>& operator-=(const vector_t<q, U> &v) {
      for (q_t i = 0; i < q; i++) {
        (*this)[i] -= (T)v[i];
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

    inline T operator*(const vector_t<q, T>& v) const {
      double prod = 0;

      for (q_t i = 0; i < q; i++) {
        prod += v[i] * (*this)[i];
      }

      return prod;
    }

    template <class U>
    inline vector_t<q, T> operator/(U a) const {
      vector_t<q, T> result;
      for (q_t i = 0; i < q; i++) {
        result[i] = (*this)[i] / a;
      }

      return result;
    }
};

template<q_t q, class T>
inline vector_t<q, T> operator*(v_t a, const vector_t<q, T>&v) {
  return v * a;
}

template<q_t q, class T>
inline vector_t<q, double> operator*(double a, const vector_t<q, T>&v) {
  return v * a;
}

template<q_t q, class T>
std::ostream& operator<<(std::ostream& os, const vector_t<q, T>&v) {
  os << "( ";
  for (T vi : v) {
    os << vi << " ";
  }
  os << ")";
  return os;
}

