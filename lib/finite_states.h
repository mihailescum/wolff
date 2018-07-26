
#pragma once

#include <cmath>
#include <functional>
#include <array>

#define FINITE_STATES

// must have N_STATES, states[N_STATES], and state_to_ind defined before
// invoking header

std::array<std::array<std::array<double, N_STATES>, N_STATES>, N_STATES> J_probs;
std::array<std::array<double, N_STATES>, N_STATES> H_probs;

template <class X_t>
void initialize_probs(std::function <double(X_t, X_t)> J, std::function <double(X_t)> H, double T) {
  for (q_t i = 0; i < N_STATES; i++) {
    for (q_t j = 0; j < N_STATES; j++) {
      for (q_t k = 0; k < N_STATES; k++) {
        J_probs[i][j][k] = 1.0 - exp(-(J(states[i], states[k]) - J(states[j], states[k])) / T);
      }
    }
  }
  for (q_t i = 0; i < N_STATES; i++) {
    for (q_t j = 0; j < N_STATES; j++) {
      H_probs[i][j] = 1.0 - exp(-(H(states[i]) - H(states[j])) / T);
    }
  }
}




