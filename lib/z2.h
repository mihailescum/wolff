
#pragma once

#include <gsl/gsl_rng.h>
#include "types.h"
#include "ising.h"
#include "state.h"

struct z2_t { bool x; };

void init(z2_t *p) {
  p->x = false;
}

void free_spin(z2_t p) {
  // do nothing!
}

z2_t copy(z2_t x) {
  return x;
}

ising_t act(z2_t r, ising_t s) {
  ising_t rs;

  if (r.x) {
    rs.x = !s.x;
    return rs;
  } else {
    rs.x = s.x;
    return rs;
  }
}

z2_t act(z2_t r1, z2_t r2) {
  z2_t r3;

  if (r1.x) {
    r3.x = !r2.x;
    return r3;
  } else {
    r3.x = r2.x;
    return r3;
  }
}

ising_t act_inverse(z2_t r, ising_t s) {
  return act(r, s);
}

z2_t act_inverse(z2_t r1, z2_t r2) {
  return act(r1, r2);
}

// these are all functions necessary for wolff.h

z2_t generate_ising_rotation(gsl_rng *r, const state_t <z2_t, ising_t> *s) {
  z2_t rot;
  rot.x = true;
  return rot;
}



