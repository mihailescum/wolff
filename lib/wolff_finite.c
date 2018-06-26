
#include "cluster_finite.h"

void wolff_finite(state_finite_t *s, count_t sweeps, count_t sweeps_per_measurement, count_t n_measurements, measurement_t *measurements) {
  for (count_t i = 0; i < sweeps; i++) {

    count_t n_flips = 0;

    while (n_flips / h->nv < sweeps_per_measurement) {
      v_t v0 = gsl_rng_uniform_int(r, h->nv);
      R_t step;
     
      bool changed = false;
      while (!changed) {
        step = gsl_rng_uniform_int(r, s->n_transformations);
        if v(symmetric_act(s->transformations + q * step, s->spins[v0]) != s->spins[v0]) {
          changed = true;
        }
      }

      v_t tmp_flips = flip_cluster_finite(s, v0, step, r);
      n_flips += tmp_flips;

      if (n_runs > 0) {
        n_steps++;
        meas_update(clust, tmp_flips);

        if (record_autocorrelation && n_steps % ac_skip == 0) {
          update_autocorr(autocorr, s->E);
        }

      }

    }

    for (q_t i = 0; i < q; i++) {
      meas_update(M[i], s->M[i]);
    }
    meas_update(E, s->E);

    q_t n_at_max = 0;
    q_t max_M_i = 0;
    v_t max_M = 0;

    for (q_t i = 0; i < q; i++) {
      if (s->M[i] > max_M) {
        n_at_max = 1;
        max_M_i = i;
        max_M = s->M[i];
      } else if (s->M[i] == max_M) {
        n_at_max++;
      }
    }

    if (record_distribution) {
      mag_dist[s->M[0]]++;
    }

    if (n_at_max == 1) {
      for (q_t i = 0; i < q; i++) {
        meas_update(sM[max_M_i][i], s->M[i]);
      }
      meas_update(sE[max_M_i], s->E);
      freqs[max_M_i]++;
    }

    diff = fabs(meas_dx(clust) / clust->x);
  }
}

