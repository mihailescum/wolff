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

template <q_t q>
class potts_t {
  public:
    q_t x;

    typedef int *M_t;
    typedef double *F_t;
};

template <q_t q>
void init(potts_t <q> *p) {
  p->x = 0;
}

template <q_t q>
void free_spin(potts_t <q> s) {
  // do nothing!
}

void free_spin(int *s) {
  free(s);
}

void free_spin(double *s) {
  free(s);
}

template <q_t q>
potts_t <q> copy(potts_t <q> s) {
  return s;
}

template <q_t q>
void add(typename potts_t<q>::M_t *s1, int a, potts_t <q> s2) {
  (*s1)[s2.x] += a;
}

template <q_t q>
void add(typename potts_t<q>::F_t *s1, double a, potts_t <q> s2) {
  (*s1)[s2.x] += a;
}

template <q_t q>
typename potts_t<q>::M_t scalar_multiple(int factor, potts_t <q> s) {
  int *M = (int *)calloc(q, sizeof(int));
  M[s.x] += factor;
  return M;
}

template <q_t q>
typename potts_t<q>::F_t scalar_multiple(double factor, potts_t <q> s) {
  double *F = (double *)calloc(q, sizeof(double));
  F[s.x] += factor;
  return F;
}

template <q_t q>
double norm_squared(typename potts_t<q>::F_t s) {
  double total = 0;
  for (q_t i = 0; i < q; i++) {
    total += pow(s[i], 2);
  }

  return total * (double)q / ((double)q - 1.0);
}

template <q_t q>
void write_magnetization(typename potts_t<q>::M_t M, FILE *outfile) {
  fwrite(&M, sizeof(int), q, outfile);
}

