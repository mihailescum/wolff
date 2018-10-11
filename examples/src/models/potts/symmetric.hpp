
#pragma once

#include <stdlib.h>
#include <array>
#include <wolff/types.h>
#include "potts.hpp"

template <q_t q>
class symmetric_t : public std::array<q_t, q> {
  public:

    symmetric_t() {
      for (q_t i = 0; i < q; i++) {
        (*this)[i] = i;
      }
    }

    potts_t<q> act(const potts_t<q> &s) const {
      return potts_t<q>((*this)[s.x]);
    }

    symmetric_t<q> act(const symmetric_t<q>& r) const {
      symmetric_t<q> r_rot;
      for (q_t i = 0; i < q; i++) {
        r_rot[i] = (*this)[r[i]];
      }

      return r_rot;
    }

    potts_t<q> act_inverse(const potts_t<q>& s) const {
      for (q_t i = 0; i < q; i++) {
        if ((*this)[i] == s.x) {
          return potts_t<q>(i);
        }
      }

      exit(EXIT_FAILURE);
    }

    symmetric_t<q> act_inverse(const symmetric_t<q>& r) const {
      symmetric_t<q> r_rot;
      for (q_t i = 0; i < q; i++) {
        r_rot[(*this)[i]] = r[i];
      }

      return r_rot;
    }
};

