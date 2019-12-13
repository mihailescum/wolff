
#include <getopt.h>
#include <iostream>
#include <chrono>

#define WOLFF_USE_FINITE_STATES

#include <wolff_models/ising.hpp>

#include "simple_measurement.hpp"

using namespace wolff;

int main(int argc, char *argv[]) {

  // set defaults
  unsigned N = (unsigned)1e4;
  unsigned D = 2;
  unsigned L = 128;
  double T = 2.26918531421;
  double H = 0.0;

  int opt;

  // take command line arguments
  while ((opt = getopt(argc, argv, "N:D:L:T:H:")) != -1) {
    switch (opt) {
    case 'N': // number of steps
      N = (unsigned)atof(optarg);
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
  graph<> G(D, L);

  // initialize the system
  wolff::system<ising_t, ising_t, graph<>> S(G, T, Z, B);

  // initialize the random number generator
  auto seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
  std::mt19937 rng(seed);

  // define function that generates self-inverse rotations
  std::function <ising_t(std::mt19937&, const wolff::system<ising_t, ising_t, graph<>>&, const graph<>::vertex&)> gen_r = gen_ising<graph<>>;

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

