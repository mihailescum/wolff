
#include <getopt.h>

#include <wolff.h>

int main(int argc, char *argv[]) {

  count_t N = (count_t)1e7;

  D_t D = 2;
  L_t L = 128;
  double T = 2.26918531421;
  double H = 0.0;

  bool silent = false;

  int opt;

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
    case 'H': // external field
      H = atof(optarg);
      break;
    case 's': // don't print anything during simulation. speeds up slightly
      silent = true;
      break;
    default:
      exit(EXIT_FAILURE);
    }
  }

  // initialize random number generator
  gsl_rng *r = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(r, rand_seed());

  state_t <z2_t, ising_t> s(D, L, T, ising_dot, std::bind(scalar_field, std::placeholders::_1, H));

  std::function <z2_t(gsl_rng *, const state_t <z2_t, ising_t> *)> gen_R = generate_ising_rotation;

  double average_M = 0;

  typedef std::function <void(const state_t <z2_t, ising_t> *)> meas_func;

  meas_func measurement = [&] (const state_t <z2_t, ising_t> *s) {
    average_M += (double)s->M / (double)N / (double)s->nv;
  };

  wolff(N, &s, gen_R, measurement, r, silent);

  printf("%" PRIcount " Ising runs completed. D = %" PRID ", L = %" PRIL ", T = %g, H = %g, <M> = %g\n", N, D, L, T, H, average_M);

  gsl_rng_free(r);

  return 0;
}

