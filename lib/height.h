
#pragma once

#include <cmath>

#include "types.h"

// object definition
template <class T>
struct height_t { T x; };

// init, copy, add, subtract, scalar_multiple, free_spin, and
// write_magnetization are necessary for the operation of wolff.h
template <class T>
void init(height_t *ptr) {
  ptr->x = (T)0;
}

template <class T>
height_t <T> copy (height_t h) {
  return h;
}

template <class T>
void add (height_t <T> *h1, height_t <T> h2) {
  h1->x += h2.x;
}

template <class T>
void subtract (height_t <T> *h1, height_T <T> h2) {
  h1->x -= h2.x;
}

template <class T>
height_t <T> scalar_multiple(v_t a, height_t <T> h) {
  height_t <T> hm;
  hm.x = a * h.x;

  return hm;
}

template <class T>
void free_spin (height_t <T> h) {
}

template <class T>
void write_magnetization(height_t <T> M, FILE *outfile) {
  fwrite(&(M.x), sizeof(T), 1, outfile);
}

template <class T>
double correlation_component(height_t <T> h) {
  return (double)h.x;
}

// below here are not necessary for operation

template <class T>
T dot(height_t <T> h1, height_t <T> h2) {
  return (h1.x - h2.x) * (h1.x - h2.x);
}

