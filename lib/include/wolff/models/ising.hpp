
#ifndef WOLFF_MODELS_ISING
#define WOLFF_MODELS_ISING

#define WOLFF_FINITE_STATES_N 2

#include "../types.h"
#include "../system.hpp"

namespace wolff {

class ising_t {
  public:
    bool x;

    ising_t() : x(false) {}

    ising_t(bool x) : x(x) {} 
    ising_t(q_t x) : x((bool)x) {}

    ising_t act(const ising_t& s) const {
      if (x) {
        return ising_t(!s.x);
      } else {
        return ising_t(s.x);
      }
    }

    ising_t act_inverse(const ising_t& s) const {
      return this->act(s);
    }

    typedef int M_t;
    typedef double F_t;

    inline M_t operator*(v_t a) const {
      if (x) {
        return -(M_t)a;
      } else {
        return (M_t)a;
      }
    }

    inline F_t operator*(double a) const {
      if (x) {
        return -a;
      } else {
        return a;
      }
    }

    inline M_t operator-(const ising_t &s) const {
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

    inline int operator*(const ising_t& s) const {
      if (x == s.x) {
        return 1;
      } else {
        return -1;
      }
    }

    q_t enumerate() const {
      return (q_t)x;
    }
};

inline ising_t::M_t operator*(v_t a, const ising_t& s) {
  return s * a;
}

inline ising_t::F_t operator*(double a, const ising_t& s) {
  return s * a;
}

ising_t gen_ising(std::mt19937&, const system<ising_t, ising_t>&, v_t) {
  return ising_t(true);
};

}

#endif

