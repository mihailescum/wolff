

#ifndef WOLFF_FINITE_STATES
#define WOLFF_FINITE_STATES

#include <cmath>
#include <array>

#include "system.hpp"

namespace wolff {

std::array<std::array<std::array<double, WOLFF_FINITE_STATES_N>, WOLFF_FINITE_STATES_N>, WOLFF_FINITE_STATES_N> finite_states_Zp;
#ifndef WOLFF_NO_FIELD
std::array<std::array<double, WOLFF_FINITE_STATES_N>, WOLFF_FINITE_STATES_N> finite_states_Bp;
#endif

template <class R_t, class X_t>
void finite_states_init(const system<R_t, X_t>& S) {
#ifndef WOLFF_NO_FIELD
  for (q_t i = 0; i < WOLFF_FINITE_STATES_N; i++) {
    for (q_t j = 0; j < WOLFF_FINITE_STATES_N; j++) {
      finite_states_Bp[i][j] = 1.0 - exp(-(S.B(finite_states_possible[i]) - S.B(finite_states_possible[j])) / S.T);
    }
  }
#endif
  for (q_t i = 0; i < WOLFF_FINITE_STATES_N; i++) {
    for (q_t j = 0; j < WOLFF_FINITE_STATES_N; j++) {
      for (q_t k = 0; k < WOLFF_FINITE_STATES_N; k++) {
        finite_states_Zp[i][j][k] = 1.0 - exp(-(S.Z(finite_states_possible[i], finite_states_possible[k]) - S.Z(finite_states_possible[j], finite_states_possible[k])) / S.T);
      }
    }
  }
}

}

#endif

