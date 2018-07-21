
#pragma once

#include <stdlib.h>
#include <cmath>

#include "types.h"

/* The following is the minimum definition of a spin class.
 *
 * The class must contain an M_t and an F_t for holding the sum of an
 * integer number of spins and a double-weighted number of spins,
 * respectively.
 *
 * void init(X_t *p);
 * void free_spin(X_t p);
 * X_t copy(X_t x);
 * void add(M_t *x1, int factor, X_t x2);
 * void add(F_t *x1, double factor, X_t x2);
 * M_t scalar_multiple(int factor, X_t x);
 * double norm_squared(F_t x);
 * void write_magnetization(M_t M, FILE *outfile);
 *
 */

template <q_t q, class T>
class vector_t {
  public: 
    T *x;

    // M_t needs to hold the sum of nv spins
    typedef vector_t <q, T> M_t;

    // F_t needs to hold the double-weighted sum of spins
    typedef vector_t <q, double> F_t;
};

template <q_t q, class T>
void init(vector_t <q, T> *ptr) {
  ptr->x = (T *)calloc(q, sizeof(T));

  // initialize a unit vector
  ptr->x[0] = (T)1;
}

template <q_t q, class T>
void free_spin (vector_t <q, T> v) {
  free(v.x);
}

template <q_t q, class T>
vector_t <q, T> copy (vector_t <q, T> v) {
  vector_t <q, T> v_copy;
 
 v_copy.x = (T *)calloc(q, sizeof(T));

  for (q_t i = 0; i < q; i++) {
    v_copy.x[i] = v.x[i];
  }

  return v_copy;
}

template <q_t q, class T, class U, class V>
void add(vector_t<q, U> *v1, V a, vector_t <q, T> v2) {
  for (q_t i = 0; i < q; i++) {
    v1->x[i] += (U)(a * v2.x[i]);
  }
}

template <q_t q, class T>
vector_t <q, T> scalar_multiple(int a, vector_t <q, T> v) {
  vector_t <q, T> multiple;
  multiple.x = (T *)malloc(q * sizeof(T));
  for (q_t i = 0; i < q; i++) {
    multiple.x[i] = a * v.x[i];
  }

  return multiple;
}

template <q_t q, class T>
double norm_squared (vector_t <q, T> v) {
  double tmp = 0;

  for (q_t i = 0; i < q; i++) {
    tmp += pow(v.x[i], 2);
  }

  return tmp;
}

template <q_t q, class T>
void write_magnetization(vector_t <q, T> M, FILE *outfile) {
  fwrite(M.x, sizeof(T), q, outfile);
}

// below functions and definitions are unnecessary for wolff.h but useful.

template <q_t q> // save some space and don't write whole doubles
void write_magnetization(vector_t <q, double> M, FILE *outfile) {
  for (q_t i = 0; i < q; i++) {
    float M_tmp = (float)M.x[i];
    fwrite(&M_tmp, sizeof(float), 1, outfile);
  }
}

template <q_t q, class T>
T dot(vector_t <q, T> v1, vector_t <q, T> v2) {
  T prod = 0;

  for (q_t i = 0; i < q; i++) {
    prod += v1.x[i] * v2.x[i];
  }

  return prod;
}

template <q_t q, class T>
double H_vector(vector_t <q, T> v1, T *H) {
  vector_t <q, T> H_vec;
  H_vec.x = H;
  return (double)(dot <q, T> (v1, H_vec));
}

char const *ON_strings[] = {"TRIVIAL", "ISING", "PLANAR", "HEISENBERG"};

