
#pragma once

#include <cmath>
#include <stdio.h>

#include "types.h"

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
};

template <class T>
void init(height_t<T> *ptr) {
  ptr->x = (T)0;
}

template <class T>
void free_spin(height_t<T> h) {
  // do nothing!
}

template <class T>
void free_spin(T h) {
  // do nothing!
}

void free_spin(double h) {
  // do nothing!
}

template <class T>
height_t<T> copy(height_t<T> h) {
  return h;
}

template <class T>
void add(T *h1, int a, height_t<T> h2) {
  (*h1) += a * h2.x;
}

template <class T>
void add(double *h1, double a, height_t<T> h2) {
  (*h1) += a * h2.x;
}

template <class T>
T scalar_multiple(int factor, height_t<T> h) {
  return factor * h.x;
}

template <class T>
double scalar_multiple(double factor, height_t<T> h) {
  return factor * h.x;
}

double norm_squared(double h) {
  return pow(h, 2);
}

template <class T>
void write_magnetization(T M, FILE *outfile) {
  fwrite(&M, sizeof(T), 1, outfile);
}

