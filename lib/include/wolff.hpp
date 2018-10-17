
#ifndef WOLFF_H
#define WOLFF_H

#include "wolff/cluster.hpp"

template <class R_t, class X_t>
void wolff(N_t N, wolff_system<R_t, X_t>& S,
           std::function <R_t(std::mt19937&, X_t)> r_gen,
           wolff_measurement<R_t, X_t>& A, std::mt19937& rng) {

  std::uniform_int_distribution<v_t> dist(0, S.nv - 1);

  for (N_t n = 0; n < N; n++) {
    v_t i0 = dist(rng);
    R_t r = r_gen(rng, S.s[i0]);

    A.pre_cluster(n, N, S, i0, r);

    wolff_cluster_flip<R_t, X_t>(S, i0, r, rng, A);

    A.post_cluster(n, N, S);
  }
}

#endif

