
#pragma once

#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <cmath>
#include <vector>
#include <stack>

#include "types.h"
#include "state.h"
#include "graph.h"

template <class R_t, class X_t>
void flip_cluster(state_t <R_t, X_t> *state, v_t v0, R_t r, gsl_rng *rand) {
  v_t nv = 0;

  std::stack<v_t> stack;
  stack.push(v0);

  std::vector<bool> marks(state->g.nv, false);

  while (!stack.empty()) {
    v_t v = stack.top();
    stack.pop();

    if (!marks[v]) {
      X_t si_old, si_new;
      R_t R_old, R_new;

      R_old = state->R;
      marks[v] = true;

      if (v == state->nv) {
        R_new = act (r, R_old);
      } else {
        si_old = state->spins[v];
        si_new = act (r, si_old);
      }

      for (v_t vn : state->g.v_adj[v]) {
        X_t sj;

        if (vn != state->nv) {
          sj = state->spins[vn];
        }

        double prob;

        bool is_ext = (v == state->nv || vn == state->nv);

        if (is_ext) {
          X_t rs_old, rs_new;
          v_t non_ghost;
          if (vn == state->nv) {
            rs_old = act_inverse (R_old, si_old);
            rs_new = act_inverse (R_old, si_new);
            non_ghost = v;
          } else {
            rs_old = act_inverse (R_old, sj);
            rs_new = act_inverse (R_new, sj);
            non_ghost = vn;
          }
          double dE = state->H(rs_old) - state->H(rs_new);
#ifdef FINITE_STATES
          prob = H_probs[N_STATES * state_to_ind(rs_old) + state_to_ind(rs_new)];
#else
          prob = 1.0 - exp(-dE / state->T);
#endif

          state->M -= rs_old;
          state->M += rs_new;

          state->E += dE;

          for (D_t i = 0; i < state->D; i++) {
            L_t x = (non_ghost / (v_t)pow(state->L, state->D - i - 1)) % state->L;

            state->ReF[i] -= rs_old * state->precomputed_cos[x];
            state->ReF[i] += rs_new * state->precomputed_cos[x];

            state->ImF[i] -= rs_old * state->precomputed_sin[x];
            state->ImF[i] += rs_new * state->precomputed_sin[x];
          }
        } else {
          double dE = state->J(si_old, sj) - state->J(si_new, sj);
#ifdef FINITE_STATES
          prob = J_probs[N_STATES * N_STATES * state_to_ind(si_old) + N_STATES * state_to_ind(si_new) + state_to_ind(sj)];
#else
          prob = 1.0 - exp(-dE / state->T);
#endif
          state->E += dE;
        }

        if (gsl_rng_uniform(rand) < prob) { // and with probability...
          stack.push(vn); // push the neighboring vertex to the stack
        }
      }

      if (v == state->nv) {
        free_spin(state->R);
        state->R = R_new;
      } else {
        state->spins[v] = si_new;
      }

      if (v != state->nv) { // count the number of non-external sites that flip
        nv++;
      }
    }
  }

  state->last_cluster_size = nv;
}

