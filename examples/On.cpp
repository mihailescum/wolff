
#include <getopt.h>
#include <iostream>
#include <chrono>

#include "simple_measurement.hpp"

#include <wolff/models/vector.hpp>
#include <wolff/models/orthogonal.hpp>
#include <wolff.hpp>

int main(int argc, char *argv[]) {

  // set defaults
  N_t N = (N_t)1e4;
  D_t D = 2;
  L_t L = 128;
  double T = 0.8;
  vector_t<WOLFF_N, double> H;
  H.fill(0.0);
  q_t Hi = 0;

  int opt;

  // take command line arguments
  while ((opt = getopt(argc, argv, "N:D:L:T:H:")) != -1) {
    switch (opt) {
    case 'N': // number of steps
      N = (N_t)atof(optarg);
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
      H[Hi] = atof(optarg);
      Hi++;
      break;
    default:
      exit(EXIT_FAILURE);
    }
  }

  // define the spin-spin coupling
  std::function <double(const vector_t<WOLFF_N, double>&, const vector_t<WOLFF_N, double>&)> Z = [] (const vector_t<WOLFF_N, double>& s1, const vector_t<WOLFF_N, double>& s2) -> double {
    return s1 * s2;
  };

  // define the spin-field coupling
  std::function <double(const vector_t<WOLFF_N, double>&)> B = [&] (const vector_t<WOLFF_N, double>& s) -> double {
    return H * s;
  };

  // initialize the system
  wolff_system<orthogonal_t<WOLFF_N, double>, vector_t<WOLFF_N, double>> S(D, L, T, Z, B);

  std::function <orthogonal_t<WOLFF_N, double>(std::mt19937&, const vector_t<WOLFF_N, double>&)> gen_R = generate_rotation_uniform<WOLFF_N>;

  // initailze the measurement object
  simple_measurement A(S);

  // initialize the random number generator
  auto seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
  std::mt19937 rng{seed};

  // run wolff N times
  wolff<orthogonal_t<WOLFF_N, double>, vector_t<WOLFF_N, double>>(N, S, gen_R, A, rng);

  // print the result of our measurements
  std::cout << "Wolff complete!\nThe average energy per site was " << A.avgE() / S.nv
    << ".\nThe average magnetization per site was " << A.avgM() / S.nv
    << ".\nThe average cluster size per site was " << A.avgC() / S.nv << ".\n";

  // exit
  return 0;
}

