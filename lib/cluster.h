
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
      X_t si_new;
      R_t R_new;

      marks[v] = true;

      bool v_is_ghost = (v == s.nv); // ghost site has the last index

      if (v_is_ghost) {
        R_new = r.act(s.R); // if we are, then we're moving the transformation
      } else {
        si_new = r.act(s.spins[v]); // otherwise, we're moving the spin at our site
      }

      for (const v_t &vn : s.g.v_adj[v]) {
        bool vn_is_ghost = (vn == s.nv); // any of our neighbors could be the ghost

        double dE, prob;

        if (v_is_ghost || vn_is_ghost) { // this is a ghost-involved bond
          X_t rs_old, rs_new;
          v_t non_ghost;

          if (vn_is_ghost) {
            // if our neighbor is the ghost, the current site is a normal
            // spin - rotate it back!
            rs_old = s.R.act_inverse(s.spins[v]);
            rs_new = s.R.act_inverse(si_new);
            non_ghost = v;
          } else {
            /* if we're the ghost, we need to rotate our neighbor back in
               both the old and new ways */
            rs_old = s.R.act_inverse(s.spins[vn]);
            rs_new = R_new.act_inverse(s.spins[vn]);
            non_ghost = vn;
          }

          dE = s.H(rs_old) - s.H(rs_new);

#ifdef FINITE_STATES
          prob = H_probs[state_to_ind(rs_old)][state_to_ind(rs_new)];
#endif

          s.update_magnetization(rs_old, rs_new);
          s.update_fourierZero(non_ghost, rs_old, rs_new);
        } else { // this is a perfectly normal bond!
          dE = s.J(s.spins[v], s.spins[vn]) - s.J(si_new, s.spins[vn]);

#ifdef FINITE_STATES
          prob = J_probs[state_to_ind(s.spins[v])][state_to_ind(si_new)][state_to_ind(s.spins[vn])];
#endif
        }

        s.update_energy(dE);

#ifndef FINITE_STATES
        prob = 1.0 - exp(-dE / s.T);
#endif

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

