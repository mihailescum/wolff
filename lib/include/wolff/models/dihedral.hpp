
#ifndef WOLFF_MODELS_DIHEDRAL_H
#define WOLFF_MODELS_DIHEDRAL_H

#include "potts.hpp"

namespace wolff {

#include "../types.h"

template <q_t q>
class dihedral_t {
  public:
    bool is_reflection;
    q_t x;

    dihedral_t() : is_reflection(false), x(0) {}
    dihedral_t(bool x, q_t y) : is_reflection(x), x(y) {}

    potts_t<q> act(const potts_t<q>& s) const {
      if (this->is_reflection) {
        return potts_t<q>(((q + this->x) - s.x) % q);
      } else {
        return potts_t<q>((this->x + s.x) % q);
      }
    }

    dihedral_t<q> act(dihedral_t<q> r) const {
      if (this->is_reflection) {
        return dihedral_t<q>(!(r.is_reflection), ((q + this->x) - r.x) % q); 
      } else {
        return dihedral_t<q>(r.is_reflection, (this->x + r.x) % q);
      }
    }

    potts_t<q> act_inverse(potts_t<q> s) const {
      if (this->is_reflection) {
        return this->act(s);
      } else {
        return potts_t<q>(((s.x + q) - this->x) % q);
      }
    }

    dihedral_t<q> act_inverse(dihedral_t<q> r) const {
      if (this->is_reflection) {
        return this->act(r);
      } else {
        return dihedral_t<q>(r.is_reflection, ((r.x + q) - this->x) % q);
      }
    }
};

}

#endif

