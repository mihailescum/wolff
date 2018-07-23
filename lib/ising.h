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

class ising_t {
  public:
    bool x;

    typedef int M_t;
    typedef double F_t;
};

void init(ising_t *p) {
  p->x = false;
}

void free_spin(ising_t s) {
  // do nothing!
}

void free_spin(int s) {
  // do nothing
}

void free_spin(double s) {
  // do nothing
}

ising_t copy(ising_t s) {
  return s;
}

void add(int *s1, int a, ising_t s2) {
  if (s2.x) {
    *s1 -= a;
  } else {
    *s1 += a;
  }
}

void add(double *s1, double a, ising_t s2) {
  if (s2.x) {
    *s1 -= a;
  } else {
    *s1 += a;
  }
}

int scalar_multiple(int factor, ising_t s) {
  if (s.x) {
    return -factor;
  } else {
    return factor;
  }
}

double scalar_multiple(double factor, ising_t s) {
  if (s.x) {
    return -factor;
  } else {
    return factor;
  }
}

double norm_squared(double s) {
  return pow(s, 2);
}

void write_magnetization(int M, FILE *outfile) {
  fwrite(&M, sizeof(int), 1, outfile);
}

