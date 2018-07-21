
#include <getopt.h>

#include <wolff.h>
#include <correlation.h>
#include <measure.h>

int main(int argc, char *argv[]) {

  count_t N = (count_t)1e7;

  D_t D = 2;
  L_t L = 128;
  double T = 2.26918531421;
  double H = 0.0;

  bool silent = false;
  bool N_is_sweeps = false;

  int opt;

  unsigned char measurement_flags = measurement_energy | measurement_clusterSize;

  while ((opt = getopt(argc, argv, "N:q:D:L:T:J:H:sM:S")) != -1) {
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
      H = atof(optarg);
      break;
    case 's': // don't print anything during simulation. speeds up slightly
      silent = true;
      break;
    case 'M':
      measurement_flags ^= 1 << atoi(optarg);
      break;
    case 'S':
      N_is_sweeps = true;
      break;
    default:
      exit(EXIT_FAILURE);
    }
  }

  unsigned long timestamp;

  {
    struct timespec spec;
    clock_gettime(CLOCK_REALTIME, &spec);
    timestamp = spec.tv_sec*1000000000LL + spec.tv_nsec;
  }

  FILE *outfile_info = fopen("wolff_metadata.txt", "a");

  fprintf(outfile_info, "<| \"ID\" -> %lu, \"MODEL\" -> \"ISING\", \"q\" -> 2, \"D\" -> %" PRID ", \"L\" -> %" PRIL ", \"NV\" -> %" PRIv ", \"NE\" -> %" PRIv ", \"T\" -> %.15f, \"H\" -> %.15f |>\n", timestamp, D, L, (v_t)pow(L, D), D * (v_t)pow(L, D), T, H);

  fclose(outfile_info);

  unsigned int n_measurements = 0;
  std::function <void(const state_t <z2_t, ising_t> *)> *measurements = (std::function <void(const state_t <z2_t, ising_t> *)> *)calloc(POSSIBLE_MEASUREMENTS, sizeof(std::function <void(const state_t <z2_t, ising_t> *)>));
  FILE *outfile_M, *outfile_E, *outfile_S, *outfile_F;

  if (measurement_flags & measurement_energy) {
    char *filename_E = (char *)malloc(255 * sizeof(char));
    sprintf(filename_E, "wolff_%lu_E.dat", timestamp);
    outfile_E = fopen(filename_E, "wb");
    free(filename_E);
    measurements[n_measurements] = measurement_energy_file<z2_t, ising_t> (outfile_E);
    n_measurements++;
  }

  if (measurement_flags & measurement_clusterSize) {
    char *filename_S = (char *)malloc(255 * sizeof(char));
    sprintf(filename_S, "wolff_%lu_S.dat", timestamp);
    outfile_S = fopen(filename_S, "wb");
    free(filename_S);
    measurements[n_measurements] = measurement_cluster_file<z2_t, ising_t> (outfile_S);
    n_measurements++;
  }

  if (measurement_flags & measurement_magnetization) {
    char *filename_M = (char *)malloc(255 * sizeof(char));
    sprintf(filename_M, "wolff_%lu_M.dat", timestamp);
    outfile_M = fopen(filename_M, "wb");
    free(filename_M);
    measurements[n_measurements] = measurement_magnetization_file<z2_t, ising_t> (outfile_M);
    n_measurements++;
  }

  if (measurement_flags & measurement_fourierZero) {
    char *filename_F = (char *)malloc(255 * sizeof(char));
    sprintf(filename_F, "wolff_%lu_F.dat", timestamp);
    outfile_F = fopen(filename_F, "wb");
    free(filename_F);
    measurements[n_measurements] = measurement_fourier_file<z2_t, ising_t> (outfile_F);
    n_measurements++;
  }

  meas_t *meas_sweeps;
  if (N_is_sweeps) {
    meas_sweeps = (meas_t *)calloc(1, sizeof(meas_t));
    measurements[n_measurements] = measurement_average_cluster<z2_t, ising_t> (meas_sweeps);
    n_measurements++;
  }

  // initialize random number generator
  gsl_rng *r = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(r, rand_seed());

  state_t <z2_t, ising_t> s(D, L, T, ising_dot, std::bind(scalar_field, std::placeholders::_1, H));

  std::function <z2_t(gsl_rng *, const state_t <z2_t, ising_t> *)> gen = generate_ising_rotation;

  if (N_is_sweeps) {
    count_t N_rounds = 0;
    printf("\n");
    while (N_rounds * N * meas_sweeps->x < N * s.nv) {
      printf("\033[F\033[J\033[F\033[JWOLFF: sweep %" PRIu64 " / %" PRIu64 ": E = %.2f, S = %" PRIv "\n", (count_t)(N_rounds * N * meas_sweeps->x / s.nv), N, s.E, s.last_cluster_size);
      wolff(N, &s, gen, n_measurements, measurements, r, silent);
      N_rounds++;
    }
    printf("\033[F\033[J\033[F\033[JWOLFF: sweep %" PRIu64 " / %" PRIu64 ": E = %.2f, S = %" PRIv "\n\n", (count_t)(N_rounds * N * meas_sweeps->x / s.nv), N, s.E, s.last_cluster_size);
  } else {
    wolff(N, &s, gen, n_measurements, measurements, r, silent);
  }

  free(measurements);

  if (measurement_flags & measurement_energy) {
    fclose(outfile_E);
  }
  if (measurement_flags & measurement_clusterSize) {
    fclose(outfile_S);
  }
  if (measurement_flags & measurement_magnetization) {
    fclose(outfile_M);
  }
  if (measurement_flags & measurement_fourierZero) {
    fclose(outfile_F);
  }

  if (N_is_sweeps) {
    free(meas_sweeps);
  }

  gsl_rng_free(r);

  return 0;
}

