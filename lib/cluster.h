
#pragma once

#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <cmath>

#include "types.h"
#include "state.h"
#include "stack.h"
#include "graph.h"

template <class R_t, class X_t>
void flip_cluster(state_t <R_t, X_t> *state, v_t v0, R_t r, gsl_rng *rand) {
  v_t nv = 0;

  ll_t *stack = NULL;     // create a new stack
  stack_push(&stack, v0); // push the initial vertex to the stack

  bool *marks = (bool *)calloc(state->nv + 1, sizeof(bool));

  while (stack != NULL) {
    v_t v = stack_pop(&stack);

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

      v_t nn = state->g->v_i[v + 1] - state->g->v_i[v];

      for (v_t i = 0; i < nn; i++) {
        v_t vn = state->g->v_adj[state->g->v_i[v] + i];

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

          add(&(state->M), -1, rs_old);
          add(&(state->M),  1, rs_new);
          state->E += dE;

          for (D_t i = 0; i < state->D; i++) {
            L_t x = (non_ghost / (v_t)pow(state->L, state->D - i - 1)) % state->L;

            add(&(state->ReF[i]), -state->precomputed_cos[i], rs_old);
            add(&(state->ReF[i]),  state->precomputed_cos[i], rs_new);

            add(&(state->ImF[i]), -state->precomputed_sin[i], rs_old);
            add(&(state->ImF[i]),  state->precomputed_sin[i], rs_new);
          }

          free_spin (rs_old);
          free_spin (rs_new);
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
          stack_push(&stack, vn); // push the neighboring vertex to the stack
        }
      }

      if (v == state->g->nv - 1) {
        free_spin(state->R);
        state->R = R_new;
      } else {
        free_spin(state->spins[v]);
        state->spins[v] = si_new;
      }

      if (v != state->g->nv - 1) { // count the number of non-external sites that flip
        nv++;
      }
    }
  }

  free(marks);

  state->last_cluster_size = nv;
}

