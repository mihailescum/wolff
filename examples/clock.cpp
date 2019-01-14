
#include <getopt.h>
#include <iostream>
#include <chrono>

#define WOLFF_USE_FINITE_STATES
#define WOLFF_FINITE_STATES_N WOLFF_POTTSQ
#include <wolff_models/potts.hpp>
#include <wolff_models/dihedral.hpp>
#include <wolff.hpp>

#include "simple_measurement.hpp"

using namespace wolff;

int main(int argc, char *argv[]) {

  // set defaults
  unsigned N = (unsigned)1e4;
  unsigned D = 2;
  unsigned L = 128;
  double T = 2.26918531421;
  vector_t<2, double> H;
  H.fill(0.0);
  unsigned Hi = 0;

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
      H[Hi] = atof(optarg);
      Hi++;
      break;
    default:
      exit(EXIT_FAILURE);
    }
  }

  // define the spin-spin coupling
  std::function <double(const potts_t<WOLFF_POTTSQ>&, const potts_t<WOLFF_POTTSQ>&)> Z = [] (const potts_t<WOLFF_POTTSQ>& s1, const potts_t<WOLFF_POTTSQ>& s2) -> double {
    return cos(2 * M_PI * (double)(s1.x + WOLFF_POTTSQ - s2.x) / (double)WOLFF_POTTSQ);
  };

  // define the spin-field coupling
  std::function <double(const potts_t<WOLFF_POTTSQ>&)> B = [=] (const potts_t<WOLFF_POTTSQ>& s) -> double {
    return H[0] * cos(2 * M_PI * (double)s.x / (double)WOLFF_POTTSQ) + H[1] * sin(2 * M_PI * (double)s.x / (double)WOLFF_POTTSQ);
  };

  // initialize the lattice
  graph<> G(D, L);

  // initialize the system
  system<dihedral_t<WOLFF_POTTSQ>, potts_t<WOLFF_POTTSQ>, graph<>> S(G, T, Z, B);

  // initialize the random number generator
  auto seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
  std::mt19937 rng{seed};

  // define function that generates self-inverse rotations
  std::function <dihedral_t<WOLFF_POTTSQ>(std::mt19937&, const system<dihedral_t<WOLFF_POTTSQ>, potts_t<WOLFF_POTTSQ>, graph<>>&, const graph<>::vertex&)> gen_r = [] (std::mt19937& r, const system<dihedral_t<WOLFF_POTTSQ>, potts_t<WOLFF_POTTSQ>, graph<>>& S, const graph<>::vertex& v) -> dihedral_t<WOLFF_POTTSQ> {
    dihedral_t<WOLFF_POTTSQ> rot;
    rot.is_reflection = true;
    std::uniform_int_distribution<unsigned> dist(0, WOLFF_POTTSQ - 2);
    unsigned x = dist(r);
    rot.x = (2 * S.s[v.ind].x + x + 1) % WOLFF_POTTSQ;

    return rot;
  };

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

