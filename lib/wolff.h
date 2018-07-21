
#include "cluster.h"
#include "state.h"

template <class R_t, class X_t>
void wolff(count_t N, state_t <R_t, X_t> *s, std::function <R_t(gsl_rng *, const state_t <R_t, X_t> *)> gen_R, unsigned int n_measurements, std::function <void(const state_t <R_t, X_t> *)> *measurements, gsl_rng *r, bool silent) {

  if (!silent) printf("\n");
  for (count_t steps = 0; steps < N; steps++) {
    if (!silent) printf("\033[F\033[JWOLFF: step %" PRIu64 " / %" PRIu64 ": E = %.2f, S = %" PRIv "\n", steps, N, s->E, s->last_cluster_size);

    v_t v0 = gsl_rng_uniform_int(r, s->nv);
    R_t step = gen_R(r, s);
    flip_cluster <R_t, X_t> (s, v0, step, r);
    free_spin(step);

    for (unsigned int i = 0; i < n_measurements; i++) {
      measurements[i](s);
    }
  }

  if (!silent) {
    printf("\033[F\033[J");
  }
  printf("WOLFF: step %" PRIu64 " / %" PRIu64 ": E = %.2f, S = %" PRIv "\n", N, N, s->E, s->last_cluster_size);

}

