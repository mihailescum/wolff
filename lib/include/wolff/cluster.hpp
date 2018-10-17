
#ifndef WOLFF_CLUSTER_H
#define WOLFF_CLUSTER_H

#include <random>
#include <cmath>
#include <vector>
#include <queue>

#include "types.h"
#include "system.hpp"
#include "graph.hpp"
#include "measurement.hpp"

template <class R_t, class X_t>
void wolff_cluster_flip(wolff_system<R_t, X_t>& S, v_t i0, const R_t& r,
                  std::mt19937& rng, wolff_measurement<R_t, X_t>& A) {
  std::uniform_real_distribution<double> dist(0.0, 1.0);

  std::queue<v_t> queue;
  queue.push(i0);

  std::vector<bool> marks(S.G.nv, false);

  while (!queue.empty()) {
    v_t i = queue.front();
    queue.pop();

    if (!marks[i]) { // don't reprocess anyone we've already visited!
      marks[i] = true;

      X_t si_new;
#ifndef WOLFF_NO_FIELD
      R_t s0_new;

      bool we_are_ghost = (i == S.nv);

      if (we_are_ghost) {
        s0_new = r.act(S.s0);
      } else
#endif
      {
        si_new = r.act(S.s[i]);
      }

      for (const v_t &j : S.G.adj[i]) {
        double dE, p;

#ifndef WOLFF_NO_FIELD
        bool neighbor_is_ghost = (j == S.nv);

        if (we_are_ghost || neighbor_is_ghost) {
          X_t s0s_old, s0s_new;
          v_t non_ghost;

          if (neighbor_is_ghost) {
            non_ghost = i;
            s0s_old = S.s0.act_inverse(S.s[i]);
            s0s_new = S.s0.act_inverse(si_new);
          } else {
            non_ghost = j;
            s0s_old = S.s0.act_inverse(S.s[j]);
            s0s_new = s0_new.act_inverse(S.s[j]);
          }

#ifdef WOLFF_SITE_DEPENDENCE
          dE = S.B(non_ghost, s0s_old) - S.B(non_ghost, s0s_new);
#else
          dE = S.B(s0s_old) - S.B(s0s_new);
#endif

#ifdef WOLFF_FINITE_STATES
          p = finite_states_Bp[finite_states_enum(s0s_old)]
                              [finite_states_enum(s0s_new)];
#endif

          // run measurement hooks for encountering a ghost bond
          A.ghost_bond_visited(S, non_ghost, s0s_old, s0s_new, dE);
        } else // this is a perfectly normal bond!
#endif
        {
#ifdef WOLFF_BOND_DEPENDENCE
          dE = S.Z(i, S.s[i], j, S.s[j]) - S.Z(i, si_new, j, S.s[j]);
#else
          dE = S.Z(S.s[i], S.s[j]) - S.Z(si_new, S.s[j]);
#endif

#ifdef WOLFF_FINITE_STATES
          p = finite_states_Zp[finite_states_enum(S.s[i])]
                              [finite_states_enum(si_new)]
                              [finite_states_enum(S.s[j])];
#endif

          // run measurement hooks for encountering a plain bond
          A.plain_bond_visited(S, i, si_new, j, dE);
        }

#ifndef FINITE_STATES
        p = 1.0 - exp(-dE / S.T);
#endif

        if (dist(rng) < p) {
          queue.push(j); // push the neighboring vertex to the queue 
        }
      }

#ifndef WOLFF_NO_FIELD
      if (we_are_ghost) {
        A.ghost_site_transformed(S, s0_new);
        S.s0 = s0_new;
      } else
#endif
      {
        A.plain_site_transformed(S, i, si_new);
        S.s[i] = si_new;
      }
    }
  }
}

#endif

