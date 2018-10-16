
#pragma once

#include <random>
#include <cmath>
#include <vector>
#include <stack>

#include "types.h"
#include "state.hpp"
#include "graph.hpp"
#include "meas.h"

template <class R_t, class X_t>
void flip_cluster(state_t<R_t, X_t>& s, v_t v0, const R_t& r, std::mt19937& rand, wolff_measurement<R_t, X_t>& m) {
  std::uniform_real_distribution<double> dist(0.0, 1.0);

  std::stack<v_t> stack;
  stack.push(v0);

  std::vector<bool> marks(s.g.nv, false);

  while (!stack.empty()) {
    v_t v = stack.top();
    stack.pop();

    if (!marks[v]) { // don't reprocess anyone we've already visited!
      marks[v] = true;

      X_t si_new;
#ifndef NOFIELD
      R_t R_new;

      bool v_is_ghost = (v == s.nv); // ghost site has the last index

      if (v_is_ghost) {
        R_new = r.act(s.R); // if we are, then we're moving the transformation
      } else
#endif
      {
        si_new = r.act(s.spins[v]); // otherwise, we're moving the spin at our site
      }

      for (const v_t &vn : s.g.v_adj[v]) {
        double dE, prob;

#ifndef NOFIELD
        bool vn_is_ghost = (vn == s.nv); // any of our neighbors could be the ghost

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

#ifdef SITE_DEPENDENCE
          dE = s.H(non_ghost, rs_old) - s.H(non_ghost, rs_new);
#else
          dE = s.H(rs_old) - s.H(rs_new);
#endif

#ifdef FINITE_STATES
          prob = H_probs[state_to_ind(rs_old)][state_to_ind(rs_new)];
#endif

          // run measurement hooks for encountering a ghost bond
          m.ghost_bond_added(non_ghost, rs_old, rs_new, dE);
        } else // this is a perfectly normal bond!
#endif
        {
#ifdef BOND_DEPENDENCE
          dE = s.J(v, s.spins[v], vn, s.spins[vn]) - s.J(v, si_new, vn, s.spins[vn]);
#else
          dE = s.J(s.spins[v], s.spins[vn]) - s.J(si_new, s.spins[vn]);
#endif


#ifdef FINITE_STATES
          prob = J_probs[state_to_ind(s.spins[v])][state_to_ind(si_new)][state_to_ind(s.spins[vn])];
#endif

          // run measurement hooks for encountering a plain bond
          m.plain_bond_added(v, s.spins[v], si_new, vn, s.spins[vn], dE);
        }

#ifndef FINITE_STATES
        prob = 1.0 - exp(-dE / s.T);
#endif

        if (dist(rand) < prob) {
          stack.push(vn); // push the neighboring vertex to the stack
        }
      }

#ifndef NOFIELD
      if (v_is_ghost) {
        m.ghost_site_transformed(s.R, R_new);
        s.R = R_new;
      } else
#endif
      {
        m.plain_site_transformed(v, s.spins[v], si_new);
        s.spins[v] = si_new;
      }
    }
  }
}

