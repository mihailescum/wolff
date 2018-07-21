#pragma once

#include <cmath>
#include <stdlib.h>

#include "types.h"

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

template <class T>
void add(T *s1, T a, ising_t s2) {
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


double norm_squared(double s) {
  return pow(s, 2);
}

template <class T>
void write_magnetization(T M, FILE *outfile) {
  fwrite(&M, sizeof(T), 1, outfile);
}

// below this line is unnecessary, but convenient

double ising_dot(ising_t s1, ising_t s2) {
  if (s1.x == s2.x) {
    return 1.0;
  } else {
    return -1.0;
  }
}

double scalar_field(ising_t s, double H) {
  if (s.x) {
    return -H;
  } else {
    return H;
  }
}

