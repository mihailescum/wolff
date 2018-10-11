#pragma once

#include <cmath>
#include <stdio.h>

#include <wolff/types.h>

class ising_t {
  public:
    bool x;

    typedef int M_t;
    typedef double F_t;

    ising_t() : x(false) {}
    ising_t(bool x) : x(x) {}
    ising_t(int x) : x((bool)x) {}

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

    inline int operator-(const ising_t &s) const {
      if (x == s.x) {
        return 0;
      } else {
        if (x) {
          return -2;
        } else {
          return 2;
        }
      }
    }
};

double norm_squared(double s) {
  return pow(s, 2);
}

void write_magnetization(int M, FILE *outfile) {
  fwrite(&M, sizeof(int), 1, outfile);
}

#define N_STATES 2
const ising_t states[2] = {ising_t(0), ising_t(1)};
q_t state_to_ind(ising_t state) { return (q_t)state.x; }

