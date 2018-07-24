
#include <stdbool.h>
#include <stdlib.h>

#include "types.h"

typedef struct {
  q_t i;
  bool r;
} dihedral_t;

dihedral_t *dihedral_compose(q_t q, q_t gti, const dihedral_t *g2);

q_t dihedral_act(q_t q, q_t gi, bool r, q_t s);

q_t dihedral_inverse_act(q_t q, const dihedral_t *g, q_t s);

q_t *dihedral_gen_transformations(q_t q);
R_t *dihedral_gen_involutions(q_t q);

R_t factorial(q_t);

#ifdef __cplusplus

template <class T, q_t q>
struct dihedral2_t { bool is_reflection; T x; };

template <class T, q_t q>
void init(dihedral2_t<T, q> *ptr) {
  ptr->is_reflection = false;
  ptr->x = (T)0;
}

template <class T, q_t q>
dihedral2_t<T, q> copy(dihedral2_t<T, q> r) {
  return r;
}

template <class T, q_t q>
void free_spin(dihedral2_t<T, q> r) {
  // do nothing!
}

template <q_t q>
potts_t<q> act(dihedral2_t<q_t, q> r, potts_t<q> s) {
  potts_t<q> s2;
  if (r.is_reflection) {
    s2.x = ((q + r.x) - s.x) % q;
  } else {
    s2.x = (r.x + s.x) % q;
  }

  return s2;
}

template <q_t q>
dihedral2_t<q_t,q> act(dihedral2_t<q_t,q> r1, dihedral2_t<q_t,q> r2) {
  dihedral2_t<q_t,q> r3;

  if (r1.is_reflection) {
    r3.is_reflection = !(r2.is_reflection);
    r3.x = ((q + r1.x) - r2.x) % q;
  } else {
    r3.is_reflection = r2.is_reflection;
    r3.x = (r1.x + r2.x) % q;
  }

  return r3;
}

template <q_t q>
potts_t<q> act_inverse(dihedral2_t<q_t,q> r, potts_t<q> s) {
  potts_t<q> s2;
  if (r.is_reflection) {
    s2.x = ((r.x + q) - s.x) % q;
  } else {
    s2.x = ((s.x + q) - r.x) % q;
  }

  return s2;
}

template <q_t q>
dihedral2_t<q_t, q> act_inverse(dihedral2_t<q_t,q> r1, dihedral2_t<q_t,q> r2) {
  dihedral2_t<q_t,q> r3;

  if (r1.is_reflection) {
    r3.is_reflection = !(r2.is_reflection);
    r3.x = ((r1.x + q) - r2.x) % q;
  } else {
    r3.is_reflection = r2.is_reflection;
    r3.x = ((r2.x + q) - r1.x) % q;
  }

  return r3;
}

#endif

