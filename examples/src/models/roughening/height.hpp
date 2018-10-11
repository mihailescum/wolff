
#pragma once

#include <cmath>
#include <stdio.h>

#include <wolff/types.h>

/* The following is the minimum definition of a spin class.
 *
 * The class must contain an M_t and an F_t for holding the sum of an
 * integer number of spins and a double-weighted number of spins,
 * respectively.
 *
 * void init(X_t *p);
 * void free_spin(X_t p);
 * void free_spin(M_t p);
 * void free_spin(F_t p);
 * X_t copy(X_t x);
 * void add(M_t *x1, int factor, X_t x2);
 * void add(F_t *x1, double factor, X_t x2);
 * M_t scalar_multiple(int factor, X_t x);
 * F_t scalar_multiple(double factor, X_t x);
 * double norm_squared(F_t x);
 * void write_magnetization(M_t M, FILE *outfile);
 *
 */

template <class T>
struct height_t {
  T x;

  typedef T M_t;
  typedef double F_t;

  height_t() : x(0) {}

  height_t(T x) : x(x) {}

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

double norm_squared(double h) {
  return pow(h, 2);
}

template <class T>
void write_magnetization(T M, FILE *outfile) {
  fwrite(&M, sizeof(T), 1, outfile);
}

