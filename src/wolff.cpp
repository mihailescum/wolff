
#include <time.h>
#include <getopt.h>

#include <cluster.h>

double H_vector(vector_t <2, double> v1, double *H) {
  vector_t <2, double> H_vec;
  H_vec.x = H;
  return dot <2, double> (v1, H_vec);
}

int main(int argc, char *argv[]) {

  count_t N = (count_t)1e7;

  D_t D = 2;
  L_t L = 128;
  double T = 2.26918531421;
  double *H = (double *)calloc(MAX_Q, sizeof(double));

  bool silent = false;

  int opt;
  q_t J_ind = 0;
  q_t H_ind = 0;

  while ((opt = getopt(argc, argv, "N:q:D:L:T:J:H:s")) != -1) {
    switch (opt) {
    case 'N': // number of steps
      N = (count_t)atof(optarg);
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

  state_t <orthogonal_t <2, double>, vector_t <2, double>> s(D, L, T, dot <2, double>, std::bind(H_vector, std::placeholders::_1, H));

  // initialize random number generator
  gsl_rng *r = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(r, rand_seed());

  unsigned long timestamp;

  {
    struct timespec spec;
    clock_gettime(CLOCK_REALTIME, &spec);
    timestamp = spec.tv_sec*1000000000LL + spec.tv_nsec;
  }

  FILE *outfile_info = fopen("wolff_metadata.txt", "a");

  fprintf(outfile_info, "<| \"ID\" -> %lu, \"D\" -> %" PRID ", \"L\" -> %" PRIL ", \"NV\" -> %" PRIv ", \"NE\" -> %" PRIv ", \"T\" -> %.15f, \"H\" -> {", timestamp, D, L, s.nv, s.ne, T);

  for (q_t i = 0; i < 2; i++) {
    fprintf(outfile_info, "%.15f", H[i]);
    if (i < 2 - 1) {
      fprintf(outfile_info, ", ");
    }
  }

  fprintf(outfile_info, "} |>\n");

  fclose(outfile_info);

  char *filename_M = (char *)malloc(255 * sizeof(char));
  char *filename_E = (char *)malloc(255 * sizeof(char));
  char *filename_S = (char *)malloc(255 * sizeof(char));

  sprintf(filename_M, "wolff_%lu_M.dat", timestamp);
  sprintf(filename_E, "wolff_%lu_E.dat", timestamp);
  sprintf(filename_S, "wolff_%lu_S.dat", timestamp);

  FILE *outfile_M = fopen(filename_M, "wb");
  FILE *outfile_E = fopen(filename_E, "wb");
  FILE *outfile_S = fopen(filename_S, "wb");

  free(filename_M);
  free(filename_E);
  free(filename_S);

  v_t cluster_size = 0;

  if (!silent) printf("\n");
  for (count_t steps = 0; steps < N; steps++) {
    if (!silent) printf("\033[F\033[JWOLFF: sweep %" PRIu64 " / %" PRIu64 ": E = %.2f, M_0 = %.2f, S = %" PRIv "\n", steps, N, s.E, s.M.x[0], cluster_size);

    v_t v0 = gsl_rng_uniform_int(r, s.nv);
     
    orthogonal_t <2, double> step;
    generate_rotation<2>(r, &step);

    cluster_size = flip_cluster <orthogonal_t <2, double>, vector_t <2, double>> (&s, v0, step, r);

    free_spin(step);

    fwrite(&(s.E), sizeof(double), 1, outfile_E);
    fwrite(s.M.x, sizeof(double), 2, outfile_M);
    fwrite(&cluster_size, sizeof(uint32_t), 1, outfile_S);

  }
  if (!silent) {
    printf("\033[F\033[J");
  }
  printf("WOLFF: sweep %" PRIu64 " / %" PRIu64 ": E = %.2f, M_0 = %.2f, S = %" PRIv "\n", N, N, s.E, s.M.x[0], cluster_size);

  fclose(outfile_M);
  fclose(outfile_E);
  fclose(outfile_S);

  gsl_rng_free(r);

  free(H);

  return 0;
}

