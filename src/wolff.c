
#include <wolff.h>

int main(int argc, char *argv[]) {
  int opt;
  bool output_state, use_dual;
  lattice_t lat;
  uint16_t L;
  uint32_t min_runs, lattice_i;
  uint64_t N;
  double T, H, eps;

  L = 128;
  N = 1000;
  lat = SQUARE_LATTICE;
  use_dual = false;
  T = 2.3;
  H = 0;
  eps = 1e30;
  output_state = false;
  min_runs = 10;

  while ((opt = getopt(argc, argv, "N:L:T:H:m:e:oq:D")) != -1) {
    switch (opt) {
    case 'N':
      N = (uint64_t)atof(optarg);
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
    case 'e':
      eps = atof(optarg);
      break;
    case 'o':
      output_state = true;
      break;
    case 'D':
      use_dual = true;
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
  s->M = h->nv;
  s->H = -(1.0 * h->ne) - H * h->nv;

  double *bond_probs = get_bond_probs(T, H, s);

  double diff = 1e31;
  uint64_t n_runs = 0;

  double E1, E2, dE1, M1, M2, dM1, C, dC, X, dX, Mmu2, Mmu4, Emu2, Emu4;
  double clust_per_sweep = 0;

  E1 = 0;
  E2 = 0;
  M1 = 0;
  M2 = 0;
  C = 0;
  dC = 0;
  X = 0;
  dX = 0;
  dE1 = 0;
  dM1 = 0;
  Mmu2 = 0;
  Mmu4 = 0;
  Emu2 = 0;
  Emu4 = 0;

  printf("\n");
  while (diff > eps || diff == 0. || n_runs < min_runs) {
    printf("\033[F\033[JWOLFF: sweep %" PRIu64
           ", dH/H = %.4f, dM/M = %.4f, dC/C = %.4f, dX/X = %.4f, cps: %.1f\n",
           n_runs, fabs(dE1 / E1), dM1 / M1, dC / C, dX / X, clust_per_sweep);

    uint32_t n_flips = 0;
    uint32_t n_clust = 0;

    while (n_flips / h->nv < N) {
      n_flips += wolff_step(T, H, s, r, bond_probs);
      n_clust++;
    }

    E1 = E1 * (n_runs / (n_runs + 1.)) + s->H * 1. / (n_runs + 1.);
    M1 = M1 * (n_runs / (n_runs + 1.)) + s->M * 1. / (n_runs + 1.);
    E2 = E2 * (n_runs / (n_runs + 1.)) + pow(s->H, 2) * 1. / (n_runs + 1.);
    M2 = M2 * (n_runs / (n_runs + 1.)) + pow(s->M, 2) * 1. / (n_runs + 1.);

    Mmu2 = Mmu2 * (n_runs / (n_runs + 1.)) +
           pow(s->M - M1, 2) * 1. / (n_runs + 1.);
    Mmu4 = Mmu4 * (n_runs / (n_runs + 1.)) +
           pow(s->M - M1, 4) * 1. / (n_runs + 1.);
    Emu2 = Emu2 * (n_runs / (n_runs + 1.)) +
           pow(s->H - E1, 2) * 1. / (n_runs + 1.);
    Emu4 = Emu4 * (n_runs / (n_runs + 1.)) +
           pow(s->H - E1, 4) * 1. / (n_runs + 1.);

    if (n_runs > 1) {
      double Msigma2 = n_runs / (n_runs - 1) * (M2 - pow(M1, 2));
      X = Msigma2 / T;
      dX =
          sqrt((Mmu4 - (n_runs - 3.) / (n_runs - 1.) * pow(Mmu2, 2)) / n_runs) /
          T;

      double Esigma2 = n_runs / (n_runs - 1) * (E2 - pow(E1, 2));
      C = Esigma2 / T;
      dC =
          sqrt((Emu4 - (n_runs - 3.) / (n_runs - 1.) * pow(Emu2, 2)) / n_runs) /
          T;

      dE1 = sqrt(Esigma2 / n_runs);
      dM1 = sqrt(Msigma2 / n_runs);

      diff = fabs(dX / X);
    }

    clust_per_sweep = clust_per_sweep * (n_runs / (n_runs + 1.)) +
                      (n_clust * 1. / N) * 1. / (n_runs + 1.);

    n_runs++;
  }
  printf("\033[F\033[JWOLFF: sweep %" PRIu64
         ", dH/H = %.4f, dM/M = %.4f, dC/C = %.4f, dX/X = %.4f, cps: %.1f\n",
         n_runs, fabs(dE1 / E1), dM1 / M1, dC / C, dX / X, clust_per_sweep);

  FILE *outfile = fopen("out.dat", "a");
  fprintf(outfile,
          "%u %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f\n", L,
          T, H, E1 / h->nv, dE1 / h->nv, M1 / h->nv, dM1 / h->nv, C / h->nv,
          dC / h->nv, X / h->nv, dX / h->nv);
  fclose(outfile);

  free(bond_probs);

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
  free(bond_probs);
  graph_free(h);

  return 0;
}
