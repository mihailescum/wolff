
#pragma once

#include <stdlib.h>
#include <cmath>

#include "types.h"

template <q_t q, class T>
struct vector_t { T *x; };

template <q_t q, class T>
void init(vector_t <q, T> *ptr) {
  ptr->x = (T *)calloc(q, sizeof(T));

  ptr->x[0] = (T)1;
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

template <q_t q, class T>
void add (vector_t <q, T> *v1, vector_t <q, T> v2) {
  for (q_t i = 0; i < q; i++) {
    v1->x[i] += v2.x[i];
  }
}

template <q_t q, class T>
void subtract (vector_t <q, T> *v1, vector_t <q, T> v2) {
  for (q_t i = 0; i < q; i++) {
    v1->x[i] -= v2.x[i];
  }
}

template <q_t q, class T>
vector_t <q, T> scalar_multiple(v_t a, vector_t <q, T> v) {
  vector_t <q, T> multiple;
  multiple.x = (T *)malloc(q * sizeof(T));
  for (q_t i = 0; i < q; i++) {
    multiple.x[i] = a * v.x[i];
  }

  return multiple;
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
void free_spin (vector_t <q, T> v) {
  free(v.x);
}

template <q_t q, class T>
void write_magnetization(vector_t <q, T> M, FILE *outfile) {
  fwrite(M.x, sizeof(T), q, outfile);
}

template <q_t q> // save some space and don't write whole doubles
void write_magnetization(vector_t <q, double> M, FILE *outfile) {
  for (q_t i = 0; i < q; i++) {
    float M_tmp = (float)M.x[i];
    fwrite(&M_tmp, sizeof(float), 1, outfile);
  }
}

template <q_t q, class T>
double correlation_component(vector_t <q, T> v) {
  return (double)v.x[0];
}

template <q_t q, class T>
double H_vector(vector_t <q, T> v1, T *H) {
  vector_t <q, T> H_vec;
  H_vec.x = H;
  return (double)(dot <q, T> (v1, H_vec));
}

char const *ON_strings[] = {"TRIVIAL", "ISING", "PLANAR", "HEISENBERG"};

