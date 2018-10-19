
#ifndef WOLFF_H
#define WOLFF_H

#include "wolff/cluster.hpp"
#include "wolff/measurement.hpp"

namespace wolff{

template <class R_t, class X_t>
void system<R_t, X_t>::run_wolff(N_t N,
           std::function <R_t(std::mt19937&, const system<R_t, X_t>&, v_t)> r_gen,
           measurement<R_t, X_t>& A, std::mt19937& rng) {

  std::uniform_int_distribution<v_t> dist(0, nv - 1);

  for (N_t n = 0; n < N; n++) {
    v_t i0 = dist(rng);
    R_t r = r_gen(rng, *this, i0);

    A.pre_cluster(n, N, *this, i0, r);

    this->flip_cluster(i0, r, rng, A);

    A.post_cluster(n, N, *this);
  }
}

}

#endif

