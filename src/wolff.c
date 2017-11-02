
#include <wolff.h>

int main(int argc, char *argv[]) {
  int opt;
  sim_t sim;
  bool output_state, use_dual, M_stop, record_autocorrelation;
  lattice_t lat;
  uint16_t L;
  uint32_t min_runs, lattice_i, sim_i;
  uint64_t N, n, W, ac_skip;
  double T, H, eps;

  sim = WOLFF;
  L = 128;
  N = 100000000000000;
  W = 10;
  ac_skip = 1;
  n = 3;
  lat = SQUARE_LATTICE;
  use_dual = false;
  M_stop = false;
  record_autocorrelation = false;
  T = 2.3;
  H = 0;
  eps = 0;
  output_state = false;
  min_runs = 10;

  while ((opt = getopt(argc, argv, "n:N:L:T:H:m:S:e:oq:DMas:W:")) != -1) {
    switch (opt) {
    case 'N':
      N = (uint64_t)atof(optarg);
      break;
    case 'W':
      W = (uint64_t)atof(optarg);
      break;
    case 'S':
      ac_skip = (uint64_t)atof(optarg);
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

  meas_t *M, *aM, *eM, *mM, *E, *eE, *mE, *clust;

  M = calloc(1, sizeof(meas_t));
  aM = calloc(1, sizeof(meas_t));
  eM = calloc(1, sizeof(meas_t));
  mM = calloc(1, sizeof(meas_t));
  E = calloc(1, sizeof(meas_t));
  eE = calloc(1, sizeof(meas_t));
  mE = calloc(1, sizeof(meas_t));
  clust = calloc(1, sizeof(meas_t));

  autocorr_t *autocorr;
  if (record_autocorrelation) {
    autocorr = (autocorr_t *)calloc(1, sizeof(autocorr_t));
    autocorr->W = W;
    autocorr->OO = (double *)calloc(W, sizeof(double));
  }

  printf("\n");
  while (((diff > eps || diff != diff) && n_runs < N) || n_runs < min_runs) {
    printf("\033[F\033[JWOLFF: sweep %" PRIu64
           ", dH/H = %.4f, dM/M = %.4f, dC/C = %.4f, dX/X = %.4f, cps: %.1f\n",
           n_runs, fabs(E->dx / E->x), M->dx / M->x, E->dc / E->c, M->dc / M->c, clust_per_sweep);

    uint32_t n_flips = 0;
    uint32_t n_clust = 0;

    while (n_flips / h->nv < n) {
      uint32_t tmp_flips = wolff_step(T, H, s, sim, r, bond_probs);
      n_flips += tmp_flips;
      n_clust++;
      n_steps++;

      if (record_autocorrelation && n_runs > 0) {
        if (n_steps % ac_skip == 0) {
          update_autocorr(autocorr, s->H);
        }
        update_meas(clust, tmp_flips);
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
      diff = fabs(eM->dc / eM->c);
    }

    clust_per_sweep = add_to_avg(clust_per_sweep, n_clust * 1. / n, n_runs);

    n_runs++;
  }

  printf("\033[F\033[JWOLFF: sweep %" PRIu64
         ", dH/H = %.4f, dM/M = %.4f, dC/C = %.4f, dX/X = %.4f, cps: %.1f\n",
         n_runs, fabs(E->dx / E->x), M->dx / M->x, E->dc / E->c, M->dc / M->c, clust_per_sweep);

  FILE *outfile = fopen("out.dat", "a");

  double tau = 0;
  double dtau = 0;
  if (record_autocorrelation) {
    double kappa = rho(autocorr, W - 1) / rho(autocorr, W - 2);
    double R = rho(autocorr, W - 1) / (1-kappa);

    double ttau = 0.5 + R;
    double X = 0.5;
    double Y = 0;
    for (uint64_t i = 0; i < W - 1; i++) {
      ttau += rho(autocorr, i);
      X += pow(rho(autocorr, i), 2);
      if (i > 0) {
        Y += 0.5 * pow(rho(autocorr, i) - rho(autocorr, i-1), 2);
      } else {
        Y += 0.5 * pow(rho(autocorr, i) - 1, 2);
      }
    }

    double dttau = sqrt(4.0 / autocorr->n * (pow(ttau, 2) * (W +(1+kappa)/(1-kappa)-0.5) + kappa * X / pow(1-kappa, 2) + pow(kappa, 2) * Y / pow(1-kappa, 4)));
    
    free(autocorr->OO);
    while (autocorr->Op != NULL) {
      stack_pop_d(&(autocorr->Op));
    }
    free(autocorr);
    
    tau = ttau * ac_skip * clust->x / h->nv;
    dtau = tau * sqrt(pow(dttau / ttau, 2) + pow(clust->dx / clust->x, 2));
  }

  fprintf(outfile,
          "%u %.15f %.15f %" PRIu64 " %" PRIu64 " %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f %" PRIu64 " %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f %" PRIu64 " %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f\n",
          L, T, H, n_runs, n_steps,
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
  free(clust);
  free(bond_probs);
  graph_free(h);

  return 0;
}
