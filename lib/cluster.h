
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
void flip_cluster(state_t<R_t, X_t>& s, v_t v0, const R_t& r, gsl_rng *rand) {
  v_t nv = 0;

  std::stack<v_t> stack;
  stack.push(v0);

  std::vector<bool> marks(s.g.nv, false);

  while (!stack.empty()) {
    v_t v = stack.top();
    stack.pop();

    if (!marks[v]) { // don't reprocess anyone we've already visited!
      X_t si_old, si_new;
      R_t R_old, R_new;

      R_old = s.R;
      marks[v] = true;

      bool v_is_ghost = (v == s.nv); // ghost site has the last index

      if (v_is_ghost) {
        R_new = r.act(R_old); // if we are, then we're moving the transformation
      } else {
        si_old = s.spins[v];
        si_new = r.act(si_old); // otherwise, we're moving the spin at our site
      }

      for (const v_t &vn : s.g.v_adj[v]) {
        X_t sj;
        bool vn_is_ghost = (vn == s.nv); // any of our neighbors could be the ghost

        if (!vn_is_ghost) {
          sj = s.spins[vn];
        }

        double prob;

        if (v_is_ghost || vn_is_ghost) { // if this is a ghost-involved bond...
          X_t rs_old, rs_new;
          v_t non_ghost;
          if (vn_is_ghost) {
            rs_old = R_old.act_inverse(si_old);
            rs_new = R_old.act_inverse(si_new);
            non_ghost = v;
          } else {
            rs_old = R_old.act_inverse(sj);
            rs_new = R_new.act_inverse(sj);
            non_ghost = vn;
          }

          double dE = s.H(rs_old) - s.H(rs_new);

#ifdef FINITE_STATES
          prob = H_probs[state_to_ind(rs_old)][state_to_ind(rs_new)];
#else
          prob = 1.0 - exp(-dE / s.T);
#endif

          s.M -= rs_old;
          s.M += rs_new;

          s.E += dE;

          for (D_t i = 0; i < s.D; i++) {
            L_t x = (non_ghost / (v_t)pow(s.L, s.D - i - 1)) % s.L;

            s.ReF[i] -= rs_old * s.precomputed_cos[x];
            s.ReF[i] += rs_new * s.precomputed_cos[x];

            s.ImF[i] -= rs_old * s.precomputed_sin[x];
            s.ImF[i] += rs_new * s.precomputed_sin[x];
          }
        } else { // otherwise, we're at a perfectly normal bond!
          double dE = s.J(si_old, sj) - s.J(si_new, sj);

#ifdef FINITE_STATES
          prob = J_probs[state_to_ind(si_old)][state_to_ind(si_new)][state_to_ind(sj)];
#else
          prob = 1.0 - exp(-dE / s.T);
#endif

          s.E += dE;
        }

        if (gsl_rng_uniform(rand) < prob) {
          stack.push(vn); // push the neighboring vertex to the stack
        }
      }

      if (v_is_ghost) {
        s.R = R_new;
      } else {
        s.spins[v] = si_new;
        nv++;
      }
    }
  }

  s.last_cluster_size = nv;
}

