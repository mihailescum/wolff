
#include <getopt.h>
#include <iostream>
#include <chrono>

#include "simple_measurement.hpp"

#include <wolff/models/ising.hpp>
#include <wolff/finite_states.hpp>

#include <wolff.hpp>

using namespace wolff;

int main(int argc, char *argv[]) {

  // set defaults
  N_t N = (N_t)1e4;
  D_t D = 2;
  L_t L = 128;
  double T = 2.26918531421;
  double H = 0.0;

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
      H = atof(optarg);
      break;
    default:
      exit(EXIT_FAILURE);
    }
  }

  // define the spin-spin coupling
  std::function <double(const ising_t&, const ising_t&)> Z = [] (const ising_t& s1, const ising_t& s2) -> double {
    return (double)(s1 * s2);
  };

  // define the spin-field coupling
  std::function <double(const ising_t&)> B = [=] (const ising_t& s) -> double {
    return H * s;
  };

  // initialize the lattice
  graph G(D, L);

  // initialize the system
  system<ising_t, ising_t> S(G, T, Z, B);

  // initialize the random number generator
  auto seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
  std::mt19937 rng{seed};

  // define function that generates self-inverse rotations
  std::function <ising_t(std::mt19937&, const system<ising_t, ising_t>&, v_t)> gen_r = gen_ising;

  // initailze the measurement object
  simple_measurement A(S);

  // run wolff N times
  S.run_wolff(N, gen_r, A, rng);

  // print the result of our measurements
  std::cout << "Wolff complete!\nThe average energy per site was " << A.avgE() / S.nv
    << ".\nThe average magnetization per site was " << A.avgM() / S.nv
    << ".\nThe average cluster size per site was " << A.avgC() / S.nv << ".\n";

  // exit
  return 0;
}

