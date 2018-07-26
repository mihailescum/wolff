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

    ising_t() {
      x = false;
    }

    ising_t(bool x) : x(x) {}

    inline int operator*(v_t a) const {
      if (x) {
        return -(int)a;
      } else {
        return (int)a;
      }
    }

    inline double operator*(double a) const {
      if (x) {
        return -a;
      } else {
        return a;
      }
    }
};

inline int& operator+=(int& M, const ising_t &s) {
  if (s.x) {
    M--;
  } else {
    M++;
  }

  return M;
}

inline int& operator-=(int& M, const ising_t &s) {
  if (s.x) {
    M++;
  } else {
    M--;
  }

  return M;
}

double norm_squared(double s) {
  return pow(s, 2);
}

void write_magnetization(int M, FILE *outfile) {
  fwrite(&M, sizeof(int), 1, outfile);
}

#define N_STATES 2
const ising_t states[2] = {{(bool)0}, {(bool)1}};
q_t state_to_ind(ising_t state) { return (q_t)state.x; }

