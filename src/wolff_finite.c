
#include <time.h>
#include <getopt.h>

#include <initial_finite.h>

int main(int argc, char *argv[]) {

  count_t N = (count_t)1e7;

  finite_model_t model = ISING;

  q_t q = 2;
  D_t D = 2;
  L_t L = 128;
  double T = 2.26918531421;
  double *J = (double *)calloc(MAX_Q, sizeof(double));
  J[0] = 1.0;
  double *H = (double *)calloc(MAX_Q, sizeof(double));

  bool silent = false;

  int opt;
  q_t J_ind = 0;
  q_t H_ind = 0;

  while ((opt = getopt(argc, argv, "N:t:q:D:L:T:J:H:s")) != -1) {
    switch (opt) {
    case 'N': // number of steps
      N = (count_t)atof(optarg);
      break;
    case 't': // type of simulation
      model = (finite_model_t)atoi(optarg);
      break;
    case 'q': // number of states, if relevant
      q = atoi(optarg);
      break;
    case 'D': // dimension
      D = atoi(optarg);
      break;
    case 'L': // linear size
      L = atoi(optarg);
      break;
    case 'T': // temperature 
      T = atof(optarg);
      break;
    case 'J': // couplings, if relevant. nth call couples states i and i + n
      J[J_ind] = atof(optarg);
      J_ind++;
      break;
    case 'H': // external field. nth call couples to state n
      H[H_ind] = atof(optarg);
      H_ind++;
      break;
    case 's': // don't print anything during simulation. speeds up slightly
      silent = true;
      break;
    default:
      exit(EXIT_FAILURE);
    }
  }

  state_finite_t *s;

  switch (model) {
    case ISING:
      s = initial_finite_prepare_ising(D, L, T, H); 
      break;
    case POTTS:
      s = initial_finite_prepare_potts(D, L, q, T, H);
      break;
    case CLOCK:
      s = initial_finite_prepare_clock(D, L, q, T, H);
      break;
    case DGM:
      s = initial_finite_prepare_dgm(D, L, q, T, H);
      break;
    default:
      printf("Not a valid model!\n");
      free(J);
      free(H);
      exit(EXIT_FAILURE);
  }

  free(J);
  free(H);

  // initialize random number generator
  gsl_rng *r = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(r, rand_seed());

  unsigned long timestamp = (unsigned long)time(NULL);

  char *filename_M = (char *)malloc(256 * sizeof(char));
  char *filename_B = (char *)malloc(256 * sizeof(char));
  char *filename_S = (char *)malloc(256 * sizeof(char));

  sprintf(filename_M, "wolff_%s_%" PRIq "_%" PRID "_%" PRIL "_%.8f", finite_model_t_strings[model], q, D, L, T);
  for (q_t i = 0; i < s->q; i++) {
    sprintf(filename_M + strlen(filename_M), "_%.8f", s->H[i]);
  }

  strcpy(filename_B, filename_M);
  strcpy(filename_S, filename_M);

  sprintf(filename_M + strlen(filename_M), "_%lu_M.dat", timestamp);
  sprintf(filename_B + strlen(filename_B), "_%lu_B.dat", timestamp);
  sprintf(filename_S + strlen(filename_S), "_%lu_S.dat", timestamp);

  FILE *outfile_M = fopen(filename_M, "wb");
  FILE *outfile_B = fopen(filename_B, "wb");
  FILE *outfile_S = fopen(filename_S, "wb");

  free(filename_M);
  free(filename_B);
  free(filename_S);

  v_t cluster_size = 0;

  if (!silent) printf("\n");
  for (count_t steps = 0; steps < N; steps++) {
    if (!silent) printf("\033[F\033[JWOLFF: sweep %" PRIu64 " / %" PRIu64 ": E = %.2f, B_0 = %" PRIv ", M_0 = %" PRIv ", S = %" PRIv "\n", steps, N, state_finite_energy(s), s->B[0], s->M[0], cluster_size);

    v_t v0 = gsl_rng_uniform_int(r, s->nv);
    R_t step;
     
    bool changed = false;
    while (!changed) {
      step = gsl_rng_uniform_int(r, s->n_involutions);
      if (symmetric_act(s->transformations + s->q * s->involutions[step], s->spins[v0]) != s->spins[v0]) {
        changed = true;
      }
    }

    cluster_size = flip_cluster_finite(s, v0, step, r);

    // v_t is never going to be bigger than 32 bits, but since it's specified
    // as a fast time many machines will actually have it be 64 bits. we cast
    // it down here to halve space.

    for (q_t i = 0; i < s->n_bond_types - 1; i++) { // if we know the occupation of all but one state we know the occupation of the last
      fwrite(&(s->B[i]), sizeof(uint32_t), 1, outfile_B);
    }

    for (q_t i = 0; i < s->q - 1; i++) { // if we know the occupation of all but one state we know the occupation of the last
      fwrite(&(s->M[i]), sizeof(uint32_t), 1, outfile_M); 
    }

    fwrite(&cluster_size, sizeof(uint32_t), 1, outfile_S);

  }
  if (!silent) {
    printf("\033[F\033[J");
  }
  printf("WOLFF: sweep %" PRIu64 " / %" PRIu64 ": E = %.2f, B_0 = %" PRIv ", M_0 = %" PRIv ", S = %" PRIv "\n", N, N, state_finite_energy(s), s->B[0], s->M[0], cluster_size);

  fclose(outfile_M);
  fclose(outfile_B);
  fclose(outfile_S);

  state_finite_free(s);
  gsl_rng_free(r);

  return 0;
}

