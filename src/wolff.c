
#include <wolff.h>

int main(int argc, char *argv[]) {
  int opt;
  sim_t sim;
  bool output_state, use_dual, M_stop, record_autocorrelation;
  lattice_t lat;
  uint16_t L;
  uint32_t min_runs, lattice_i, sim_i;
  uint64_t N, n;
  double T, H, eps;

  sim = WOLFF;
  L = 128;
  N = 100000000000000;
  n = 3;
  lat = SQUARE_LATTICE;
  use_dual = false;
  M_stop = false;
  record_autocorrelation = false;
  T = 2.3;
  H = 0;
  eps = 1e-30;
  output_state = false;
  min_runs = 10;

  while ((opt = getopt(argc, argv, "n:N:L:T:H:m:e:oq:DMas:")) != -1) {
    switch (opt) {
    case 'N':
      N = (uint64_t)atof(optarg);
      break;
    case 'n':
      n = (uint64_t)atof(optarg);
      break;
    case 'L':
      L = atoi(optarg);
      break;
    case 'T':
      T = atof(optarg);
      break;
    case 'H':
      H = atof(optarg);
      break;
    case 'm':
      min_runs = atoi(optarg);
      break;
    case 'M':
      M_stop = true;
      break;
    case 'e':
      eps = atof(optarg);
      break;
    case 'o':
      output_state = true;
      break;
    case 'D':
      use_dual = true;
      break;
    case 'a':
      record_autocorrelation = true;
      break;
    case 's':
      sim_i = atoi(optarg);
      switch (sim_i) {
      case 0:
        sim = WOLFF;
        break;
      case 1:
        sim = WOLFF_GHOST;
        break;
      case 2:
        sim = METROPOLIS;
        break;
      default:
        printf("lattice specifier must be 0 (VORONOI_LATTICE), 1 "
               "(DIAGONAL_LATTICE), or 2 (VORONOI_HYPERUNIFORM_LATTICE).\n");
        exit(EXIT_FAILURE);
      }
      break;
    case 'q':
      lattice_i = atoi(optarg);
      switch (lattice_i) {
      case 0:
        lat = SQUARE_LATTICE;
        break;
      case 1:
        lat = DIAGONAL_LATTICE;
        break;
      case 2:
        lat = TRIANGULAR_LATTICE;
        break;
      case 3:
        lat = VORONOI_HYPERUNIFORM_LATTICE;
        break;
      case 4:
        lat = VORONOI_LATTICE;
        break;
      default:
        printf("lattice specifier must be 0 (VORONOI_LATTICE), 1 "
               "(DIAGONAL_LATTICE), or 2 (VORONOI_HYPERUNIFORM_LATTICE).\n");
        exit(EXIT_FAILURE);
      }
      break;
    default:
      exit(EXIT_FAILURE);
    }
  }

  gsl_rng *r = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(r, jst_rand_seed());

  graph_t *h = graph_create(lat, TORUS_BOUND, L, use_dual);

  ising_state_t *s = (ising_state_t *)calloc(1, sizeof(ising_state_t));

  s->g = graph_add_ext(h);
  s->spins = (bool *)calloc(h->nv + 1, sizeof(bool));
  s->M = sign(H) * h->nv;
  s->H = -(1.0 * h->ne) - sign (H) * H * h->nv;

  double *bond_probs = (double *)malloc(2 * sizeof(double));
  bond_probs[0] = 1 - exp(-2 / T);
  bond_probs[1] = 1 - exp(-2 * fabs(H) / T);

  double diff = 1e31;
  uint64_t n_runs = 0;
  uint64_t n_steps = 0;
  double clust_per_sweep = 0;

  meas_t *M, *aM, *eM, *mM, *E, *eE, *mE, *corrE, *corrmE;

  M = calloc(1, sizeof(meas_t));
  aM = calloc(1, sizeof(meas_t));
  eM = calloc(1, sizeof(meas_t));
  mM = calloc(1, sizeof(meas_t));
  E = calloc(1, sizeof(meas_t));
  eE = calloc(1, sizeof(meas_t));
  mE = calloc(1, sizeof(meas_t));

  corrE = (meas_t *)calloc(1, sizeof(meas_t));
  corrmE = (meas_t *)calloc(1, sizeof(meas_t));
  double batch_mean_energy;
  uint64_t batch_size = exp(0.667 * log(N));
  uint64_t n_batch = 0;
  uint64_t batch_flips = 0;

  printf("\n");
  while ((diff > eps && n_steps < N) || n_runs < min_runs) {
    printf("\033[F\033[JWOLFF: sweep %" PRIu64
           ", dH/H = %.4f, dM/M = %.4f, dC/C = %.4f, dX/X = %.4f, cps: %.1f\n",
           n_runs, fabs(E->dx / E->x), M->dx / M->x, E->dc / E->c, M->dc / M->c, clust_per_sweep);

    uint32_t n_flips = 0;
    uint32_t n_clust = 0;

    while (n_flips / h->nv < n) {
      uint32_t tmp_flips = wolff_step(T, H, s, sim, r, bond_probs);
      n_flips += tmp_flips;
      n_clust++;

      n_steps += tmp_flips;

      if (n_runs > 0) {
        update_meas(corrE, s->H);
        if (batch_flips <= batch_size && batch_flips + tmp_flips > batch_size) {
          update_meas(corrmE, batch_mean_energy / n_batch);
          batch_mean_energy = 0;
          n_batch = 0;
          batch_flips = (int64_t)(batch_flips + tmp_flips) - (int64_t)batch_size;
        }
        batch_mean_energy += s->H;
        n_batch++;
        batch_flips += tmp_flips;
      }
    }

    double HH = 1;
    if (H < 0) {
      HH = -1;
    }

    update_meas(M, (double)(s->M));
    update_meas(aM, HH * fabs((double)(s->M)));
    update_meas(E, s->H);

    if (s->M > 0) {
      update_meas(eM, HH * fabs((double)(s->M)));
      update_meas(eE, s->H);
    } else {
      update_meas(mM, - HH * fabs((double)(s->M)));
      update_meas(mE, s->H);
    }

    if (M_stop) {
      diff = fabs(eM->dx / eM->x);
    } else {
      diff = fabs(M->dc / M->c);
    }

    clust_per_sweep = add_to_avg(clust_per_sweep, n_clust * 1. / n, n_runs);

    n_runs++;
  }

  printf("\033[F\033[JWOLFF: sweep %" PRIu64
         ", dH/H = %.4f, dM/M = %.4f, dC/C = %.4f, dX/X = %.4f, cps: %.1f\n",
         n_runs, fabs(E->dx / E->x), M->dx / M->x, E->dc / E->c, M->dc / M->c, clust_per_sweep);

  FILE *outfile = fopen("out.dat", "a");

  double tau = batch_size * corrmE->c / corrE->c;
  double dtau = tau * sqrt(pow(corrmE->dc / corrmE->c, 2) + pow(corrE->dc / corrE->c, 2));

  fprintf(outfile,
          "%u %.15f %.15f %" PRIu64 " %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f %" PRIu64 " %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f %" PRIu64 " %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f\n", L,
          T, H, n_runs,
          E->x / h->nv, E->dx / h->nv, M->x / h->nv, M->dx / h->nv, E->c / h->nv, E->dc / h->nv, M->c / h->nv, M->dc / h->nv,
          eE->n, eE->x / h->nv, eE->dx / h->nv, eM->x / h->nv, eM->dx / h->nv, eE->c / h->nv, eE->dc / h->nv, eM->c / h->nv, eM->dc / h->nv,
          mE->n, mE->x / h->nv, mE->dx / h->nv, mM->x / h->nv, mM->dx / h->nv, mE->c / h->nv, mE->dc / h->nv, mM->c / h->nv, mM->dc / h->nv,
          aM->x / h->nv, aM->dx / h->nv, aM->c / h->nv, aM->dc / h->nv, tau, dtau);
  fclose(outfile);

  if (output_state) {
    FILE *state_file = fopen("state.dat", "a");
    for (uint32_t i = 0; i < h->nv; i++) {
      fprintf(state_file, "%d ", s->spins[i]);
    }
    fprintf(state_file, "\n");
    fclose(state_file);
  }

  gsl_rng_free(r);
  graph_free(s->g);
  free(s->spins);
  free(s);
  free(M);
  free(aM);
  free(eM);
  free(mM);
  free(E);
  free(eE);
  free(mE);
  free(corrE);
  free(corrmE);
  free(bond_probs);
  graph_free(h);

  return 0;
}
