
#include "wolff/cluster.hpp"
#include "wolff/state.hpp"

template <class R_t, class X_t>
void wolff(count_t N, state_t <R_t, X_t>& s, std::function <R_t(std::mt19937&, X_t)> gen_R, wolff_measurement<R_t, X_t>& m, std::mt19937& r) {

#ifdef FINITE_STATES
#ifdef NOFIELD
  initialize_probs(s.J, s.T);
#else
  initialize_probs(s.J, s.H, s.T);
#endif
#endif

  std::uniform_int_distribution<v_t> dist(0, s.nv);

  for (count_t steps = 0; steps < N; steps++) {
    v_t v0 = dist(r);
    R_t step = gen_R(r, s.spins[v0]);

    m.pre_cluster(s, steps, N, v0, step);

    flip_cluster<R_t, X_t>(s, v0, step, r, m);

    m.post_cluster(s, steps, N);
  }

}

