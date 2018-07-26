
#pragma once

#include "types.h"
#include "ising.h"

/* The minimum definition for a group type R_t to act on a spin type X_t is
 * given by the following.
 *
 * void init(R_t *p);
 * void free_spin(R_t r);
 * R_t copy(R_t r);
 * X_t act(R_t r, X_t x);
 * R_t act(R_t r, R_t r);
 * X_t act_inverse(R_t r, X_t x);
 * R_t act_inverse(R_t r, R_t r);
 *
 */

class z2_t {
  public:
  bool x;

  z2_t() : x(false) {}

  z2_t(bool x) : x(x) {}

  ising_t act(const ising_t& s) {
    if (x) {
      return ising_t(!s.x);
    } else {
      return ising_t(s.x);
    }
  }

  z2_t act(const z2_t& r) {
    if (x) {
      return z2_t(!r.x);
    } else {
      return z2_t(r.x);
    }
  }

  ising_t act_inverse(const ising_t& s) {
    return this->act(s);
  }

  z2_t act_inverse(const z2_t& r) {
    return this->act(r);
  }
};


