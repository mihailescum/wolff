
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

  // define spin-spin coupling
  std::function <double(ising_t, ising_t)> Z = [] (ising_t s1, ising_t s2) -> double {
    if (s1.x == s2.x) {
      return 1.0;
    } else {
      return -1.0;
    }
  };

  // define spin-field coupling
  std::function <double(ising_t)> B = [=] (ising_t s) -> double {
    if (s.x) {
      return -H;
    } else {
      return H;
    }
  };

  // initialize state object
  state_t <z2_t, ising_t> s(D, L, T, Z, B);

  // define function that generates self-inverse rotations
  std::function <z2_t(gsl_rng *, const state_t <z2_t, ising_t> *)> gen_R = [] (gsl_rng *, const state_t <z2_t, ising_t> *) -> z2_t {
    z2_t rot;
    rot.x = true;
    return rot;
  };

  // define function that updates any number of measurements
  double average_M = 0;
  std::function <void(const state_t <z2_t, ising_t> *)> measurement = [&] (const state_t <z2_t, ising_t> *s) {
    average_M += (double)s->M / (double)N / (double)s->nv;
  };

  // run wolff for N cluster flips
  wolff(N, &s, gen_R, measurement, r, silent);

  // tell us what we found!
  printf("%" PRIcount " Ising runs completed. D = %" PRID ", L = %" PRIL ", T = %g, H = %g, <M> = %g\n", N, D, L, T, H, average_M);

  // free the random number generator
  gsl_rng_free(r);

  return 0;

}

