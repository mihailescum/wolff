
#include <wolff/types.h>
#include <cmath>
#include "height.hpp"

template <class T>
class dihedral_inf_t {
  public:
  bool is_reflection;
  T x;

  dihedral_inf_t() : is_reflection(false), x(0) {}
  dihedral_inf_t(bool x, T y) : is_reflection(x), x(y) {}

  height_t<T> act(const height_t<T>& h) const {
    if (this->is_reflection) {
      return height_t(this->x - h.x);
    } else {
      return height_t(this->x + h.x);
    }
  }

  dihedral_inf_t<T> act(const dihedral_inf_t<T>& r) const {
    if (this->is_reflection) {
      return dihedral_inf_t<T>(!r.is_reflection, this->x - r.x);
    } else {
      return dihedral_inf_t<T>(r.is_reflection, this->x + r.x);
    }
  }

  height_t<T> act_inverse(const height_t<T>& h) const {
    if (this->is_reflection) {
      return this->act(h);
    } else {
      return height_t(h.x - this->x);
    }
  }

  dihedral_inf_t<T> act_inverse(const dihedral_inf_t<T>& r) const {
    if (this->is_reflection) {
      return this->act(r);
    } else {
      return dihedral_inf_t<T>(r.is_reflection, r.x - this->x);
    }
  }
};

