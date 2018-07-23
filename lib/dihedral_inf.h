
#include "types.h"
#include <cmath>
#include "height.h"

template <class T>
struct dihedral_inf_t { bool is_reflection; T x; };

template <class T>
void init(dihedral_inf_t<T> *ptr) {
  ptr->is_reflection = false;
  ptr->x = (T)0;
}

template <class T>
dihedral_inf_t<T> copy(dihedral_inf_t<T> r) {
  return r;
}

template <class T>
void free_spin(dihedral_inf_t<T> r) {
  // do nothing!
}

template <class T>
height_t<T> act(dihedral_inf_t<T> r, height_t<T> h) {
  height_t<T> h2;
  if (r.is_reflection) {
    h2.x = r.x - h.x;
  } else {
    h2.x = r.x + h.x;
  }

  return h2;
}

template <class T>
dihedral_inf_t<T> act(dihedral_inf_t<T> r1, dihedral_inf_t<T> r2) {
  dihedral_inf_t<T> r3;

  if (r1.is_reflection) {
    r3.is_reflection = !(r2.is_reflection);
    r3.x = r1.x - r2.x;
  } else {
    r3.is_reflection = r2.is_reflection;
    r3.x = r1.x + r2.x;
  }

  return r3;
}

template <class T>
height_t<T> act_inverse(dihedral_inf_t<T> r, height_t<T> h) {
  height_t<T> h2;
  if (r.is_reflection) {
    h2.x = r.x - h.x;
  } else {
    h2.x = h.x - r.x;
  }

  return h2;
}

template <class T>
dihedral_inf_t<T> act_inverse(dihedral_inf_t<T> r1, dihedral_inf_t<T> r2) {
  dihedral_inf_t<T> r3;

  if (r1.is_reflection) {
    r3.is_reflection = !(r2.is_reflection);
    r3.x = r1.x - r2.x;
  } else {
    r3.is_reflection = r2.is_reflection;
    r3.x = r2.x - r1.x;
  }

  return r3;
}

