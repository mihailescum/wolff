
#ifndef WOLFF_MODELS_ISING
#define WOLFF_MODELS_ISING

#include "../types.h"

class ising_t {
  public:
    bool x;

    ising_t() : x(false) {}

    ising_t(bool x) : x(x) {} 
    ising_t(int x) : x((bool)x) {}

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
};

inline ising_t::M_t operator*(v_t a, const ising_t& s) {
  return s * a;
}

inline ising_t::F_t operator*(double a, const ising_t& s) {
  return s * a;
}

#define WOLFF_FINITE_STATES_N 2
const ising_t finite_states_possible[2] = {ising_t(0), ising_t(1)};
q_t finite_states_enum(ising_t state) { return (q_t)state.x; }

#endif

