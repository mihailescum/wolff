
#pragma once

#include <cmath>
#include <wolff/types.h>

namespace wolff {

template <class T>
struct height_t {
  T x;

  height_t() : x(0) {}
  height_t(T x) : x(x) {}

  typedef T M_t;
  typedef double F_t;

  inline T operator*(v_t a) const {
    return x * a;
  }

  inline double operator*(double a) const {
    return x * a;
  }

  inline T operator-(const height_t& h) const {
    return x - h.x;
  }
};

template <class T>
inline typename height_t<T>::M_t operator*(v_t a, height_t<T> h) {
  return h * a;
}

template <class T>
inline typename height_t<T>::F_t operator*(double a, height_t<T> h) {
  return h * a;
}

template <class T>
inline T& operator+=(T& M, const height_t<T> &h) {
  M += h.x;

  return M;
}

template <class T>
inline T& operator-=(T& M, const height_t<T> &h) {
  M -= h.x;

  return M;
}

}

