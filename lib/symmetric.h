
#pragma once

#include <stdlib.h>

#include "types.h"

#ifdef __cplusplus
#include "potts.h"
extern "C" {
#endif

q_t *symmetric_compose(q_t q, const q_t *g1, const q_t *g2);

q_t symmetric_act(const q_t *g, q_t s);

q_t *symmetric_invert(q_t q, const q_t *g);

q_t *symmetric_gen_transformations(q_t q);

#ifdef __cplusplus
}
#endif

#ifdef __cplusplus
template <q_t q>
class symmetric_t {
  public:
    q_t *perm;
};

template <q_t q>
void init(symmetric_t<q> *p) {
  p->perm = (q_t *)malloc(q * sizeof(q_t));

  for (q_t i = 0; i < q; i++) {
    p->perm[i] = i;
  }
}

template <q_t q>
void free_spin(symmetric_t<q> p) {
  free(p.perm);
}

template <q_t q>
symmetric_t<q> copy(symmetric_t<q> x) {
  symmetric_t<q> x2;
  x2.perm = (q_t *)malloc(q * sizeof(q_t));

  for (q_t i = 0; i < q; i++) {
    x2.perm[i] = x.perm[i];
  }

  return x2;
}

template <q_t q>
potts_t<q> act(symmetric_t<q> r, potts_t<q> s) {
  potts_t<q> s2;
  s2.x = r.perm[s.x];
  return s2;
}

template <q_t q>
symmetric_t<q> act(symmetric_t<q> r1, symmetric_t<q> r2) {
  symmetric_t<q> r3;
  r3.perm = (q_t *)malloc(q * sizeof(q_t));
  for (q_t i = 0; i < q; i++) {
    r3.perm[i] = r1.perm[r2.perm[i]];
  }

  return r3;
}

template <q_t q>
potts_t<q> act_inverse(symmetric_t<q> r, potts_t<q> s) {
  potts_t<q> s2;

  q_t i;

  for (i = 0; i < q; i++) {
    if (r.perm[i] == s.x) {
      break;
    }
  }

  s2.x = i;

  return s2;
}

template <q_t q>
symmetric_t<q> act_inverse(symmetric_t<q> r1, symmetric_t<q> r2) {
  symmetric_t<q> r3;
  r3.perm = (q_t *)malloc(q * sizeof(q_t));
  for (q_t i = 0; i < q; i++) {
    r3.perm[r1.perm[i]] = r2.perm[i];
  }

  return r3;
}
#endif

