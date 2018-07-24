
#pragma once

#include <cmath>
#include <functional>

#define FINITE_STATES

// must have N_STATES, states[N_STATES], and state_to_ind defined before
// invoking header

double J_probs[N_STATES * N_STATES * N_STATES];
double H_probs[N_STATES * N_STATES];

template <class X_t>
void initialize_probs(std::function <double(X_t, X_t)> J, std::function <double(X_t)> H, double T) {
  for (q_t i = 0; i < N_STATES; i++) {
    for (q_t j = 0; j < N_STATES; j++) {
      for (q_t k = 0; k < N_STATES; k++) {
        J_probs[i * N_STATES * N_STATES + j * N_STATES +k] = 1.0 - exp(-(J(states[i], states[k]) - J(states[j], states[k])) / T);
      }
    }
  }
  for (q_t i = 0; i < N_STATES; i++) {
    for (q_t j = 0; j < N_STATES; j++) {
      H_probs[i * N_STATES + j] = 1.0 - exp(-(H(states[i]) - H(states[j])) / T);
    }
  }
}




