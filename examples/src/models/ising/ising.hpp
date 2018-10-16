#pragma once

#include <cmath>
#include <stdio.h>

#include <wolff/types.h>

// all that is required to use wolff.hpp is a default constructor
class ising_t {
  public:
    bool x;

    ising_t() : x(false) {}

    // optional constructors for syntactic sugar
    ising_t(bool x) : x(x) {} 
    ising_t(int x) : x((bool)x) {}

    /* below this comment is code required only for using measure.hpp in the
     * examples folder, which provides an interface for measuring several
     * generic features of models. these require
     *
     *  - an M_t, representing the magnetization or sum of all spins
     *  - an F_t, representing a double-weighted version of the magnetization
     *  - the overloaded operator *, which takes a v_t (unsigned int) and returns an M_t
     *  - the overloaded operator *, which takes a double and returns an F_t
     *  - the overloaded operator -, which takes another X_t and returns an M_t
     */

    typedef int M_t;
    typedef double F_t;

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

/* using measure.hpp additionally requires a norm_squared function which takes
 * an F_t to a double, and a write_magnetization function, which takes an M_t
 * and a FILE pointer and appropriately records the contents of the former to
 * the latter.
 */

double norm_squared(double s) {
  return pow(s, 2);
}

void write_magnetization(int M, FILE *outfile) {
  fwrite(&M, sizeof(int), 1, outfile);
}

/* these definitions allow wolff/finite_states.hpp to be invoked and provide
 * much faster performance for models whose number of possible spin
 * configurations is finite.
 */

#define N_STATES 2
const ising_t states[2] = {ising_t(0), ising_t(1)};
q_t state_to_ind(ising_t state) { return (q_t)state.x; }

