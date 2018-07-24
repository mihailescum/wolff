
#pragma once

#include "types.h"
#include "potts.h"

template <class T, q_t q>
struct dihedral_t { bool is_reflection; T x; };

template <class T, q_t q>
void init(dihedral_t<T, q> *ptr) {
  ptr->is_reflection = false;
  ptr->x = (T)0;
}

template <class T, q_t q>
dihedral_t<T, q> copy(dihedral_t<T, q> r) {
  return r;
}

template <class T, q_t q>
void free_spin(dihedral_t<T, q> r) {
  // do nothing!
}

template <q_t q>
potts_t<q> act(dihedral_t<q_t, q> r, potts_t<q> s) {
  potts_t<q> s2;
  if (r.is_reflection) {
    s2.x = ((q + r.x) - s.x) % q;
  } else {
    s2.x = (r.x + s.x) % q;
  }

  return s2;
}

template <q_t q>
dihedral_t<q_t,q> act(dihedral_t<q_t,q> r1, dihedral_t<q_t,q> r2) {
  dihedral_t<q_t,q> r3;

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
potts_t<q> act_inverse(dihedral_t<q_t,q> r, potts_t<q> s) {
  potts_t<q> s2;
  if (r.is_reflection) {
    s2.x = ((r.x + q) - s.x) % q;
  } else {
    s2.x = ((s.x + q) - r.x) % q;
  }

  return s2;
}

template <q_t q>
dihedral_t<q_t, q> act_inverse(dihedral_t<q_t,q> r1, dihedral_t<q_t,q> r2) {
  dihedral_t<q_t,q> r3;

  if (r1.is_reflection) {
    r3.is_reflection = !(r2.is_reflection);
    r3.x = ((r1.x + q) - r2.x) % q;
  } else {
    r3.is_reflection = r2.is_reflection;
    r3.x = ((r2.x + q) - r1.x) % q;
  }

  return r3;
}

