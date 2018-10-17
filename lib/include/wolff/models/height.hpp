
#pragma once

#include <cmath>
#include <wolff/types.h>

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
inline T& operator+=(T& M, const height_t<T> &h) {
  M += h.x;

  return M;
}

template <class T>
inline T& operator-=(T& M, const height_t<T> &h) {
  M -= h.x;

  return M;
}

