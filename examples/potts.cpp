
#include <getopt.h>
#include <iostream>
#include <chrono>

#define WOLFF_USE_FINITE_STATES
#define WOLFF_FINITE_STATES_N WOLFF_POTTSQ

#include <wolff/models/potts.hpp>
#include <wolff/models/symmetric.hpp>

#include "simple_measurement.hpp"

#include <wolff.hpp>

using namespace wolff;

int main(int argc, char *argv[]) {

  // set defaults
  N_t N = (N_t)1e4;
  D_t D = 2;
  L_t L = 128;
  double T = 2.26918531421;
  vector_t<WOLFF_POTTSQ, double> H;
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
  std::function <double(const potts_t<WOLFF_POTTSQ>&, const potts_t<WOLFF_POTTSQ>&)> Z = [] (const potts_t<WOLFF_POTTSQ>& s1, const potts_t<WOLFF_POTTSQ>& s2) -> double {
    if (s1.x == s2.x) {
      return 1.0;
    } else {
      return 0.0;
    }
  };

  // define the spin-field coupling
  std::function <double(const potts_t<WOLFF_POTTSQ>&)> B = [=] (const potts_t<WOLFF_POTTSQ>& s) -> double {
    return H[s.x];
  };

  // initialize the lattice
  graph G(D, L);

  // initialize the system
  system<symmetric_t<WOLFF_POTTSQ>, potts_t<WOLFF_POTTSQ>> S(G, T, Z, B);

  // initialize the random number generator
  auto seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
  std::mt19937 rng{seed};

  // define function that generates self-inverse rotations
  std::function <symmetric_t<WOLFF_POTTSQ>(std::mt19937&, const system<symmetric_t<WOLFF_POTTSQ>, potts_t<WOLFF_POTTSQ>>&, v_t)> gen_r = [] (std::mt19937& r, const system<symmetric_t<WOLFF_POTTSQ>, potts_t<WOLFF_POTTSQ>>& S, v_t i0) -> symmetric_t<WOLFF_POTTSQ> {
    symmetric_t<WOLFF_POTTSQ> rot;

    std::uniform_int_distribution<q_t> dist(0, WOLFF_POTTSQ - 2);
    q_t j = dist(r);
    q_t swap_v;
    if (j < S.s[i0].x) {
      swap_v = j;
    } else {
      swap_v = j + 1;
    }

    rot[S.s[i0].x] = swap_v;
    rot[swap_v] = S.s[i0].x;

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

