
#pragma once

#include <functional>
#include <assert.h>
#include <fftw3.h>
#include <float.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <inttypes.h>
#include <cmath>
#include <stdbool.h>
#include <string.h>
#include <sys/types.h>

#include "state.h"
#include "types.h"
#include "rand.h"
#include "stack.h"
#include "convex.h"
#include "graph.h"
#include "tree.h"
#include "measurement.h"
#include "vector.h"
#include "orthogonal.h"
#include "dihedral.h"
#include "dihinf.h"
#include "yule_walker.h"

template <class T>
void init(T*);

template <class T>
T scalar_multiple(v_t a, T b);

template <class R_t, class X_t>
X_t act(R_t a, X_t b);

template <class R_t, class X_t>
X_t act_inverse(R_t a, X_t b);

template <class T>
T copy(T a);

template <class T>
void free_spin(T a);

template <class T>
T add(T, T);

template <class T>
T subtract(T, T);

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
          prob = 1.0 - exp(-dE / state->T);

          subtract (&(state->M), rs_old);
          add (&(state->M), rs_new);
          state->E += dE;

          for (D_t i = 0; i < state->D; i++) {
            double x = (double)((non_ghost / (v_t)pow(state->L, state->D - i - 1)) % state->L) / (double)state->L; 
            scalar_subtract(&(state->ReF[i]), cos(2 * M_PI * x), rs_old);
            scalar_add(&(state->ReF[i]), cos(2 * M_PI * x), rs_new);

            scalar_subtract(&(state->ImF[i]), sin(2 * M_PI * x), rs_old);
            scalar_add(&(state->ImF[i]), sin(2 * M_PI * x), rs_new);
          }

          free_spin (rs_old);
          free_spin (rs_new);
        } else {
          double dE = state->J(si_old, sj) - state->J(si_new, sj);
          prob = 1.0 - exp(-dE / state->T);
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

