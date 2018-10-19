
#include <getopt.h>
#include <iostream>
#include <chrono>

#include "simple_measurement.hpp"

#include <wolff/models/height.hpp>
#include <wolff/models/dihedral_inf.hpp>

#include <wolff.hpp>

int main(int argc, char *argv[]) {

  // set defaults
  N_t N = (N_t)1e4;
  D_t D = 2;
  L_t L = 128;
  double T = 0.8;
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
  std::function <double(const height_t<int64_t>&, const height_t<int64_t>&)> Z = [] (const height_t<int64_t>& s1, const height_t<int64_t>& s2) -> double {
    return - pow(s1.x - s2.x, 2);
  };

  // define the spin-field coupling
  std::function <double(const height_t<int64_t>&)> B = [&] (const height_t<int64_t>& s) -> double {
    return - H * pow(s.x, 2);
  };

  // initialize the lattice
  graph G(D, L);

  // initialize the system
  system<dihedral_inf_t<int64_t>, height_t<int64_t>> S(G, T, Z, B);

  bool odd_run = false;

  std::function <dihedral_inf_t<int64_t>(std::mt19937&, const system<dihedral_inf_t<int64_t>, height_t<int64_t>>&, v_t)> gen_R_IH = [&](std::mt19937& r, const system<dihedral_inf_t<int64_t>, height_t<int64_t>>& S, v_t i0) -> dihedral_inf_t<int64_t> {
    dihedral_inf_t<int64_t> rot;
    rot.is_reflection = true;

    if (odd_run) {
      std::uniform_int_distribution<v_t> dist(0, S.nv - 2);
      v_t j = i0;

      //while (S.s[j].x == S.s[i0].x) {
        v_t tmp = dist(r);

        if (tmp < i0) {
          j = tmp;
        } else {
          j = tmp + 1;
        }
      //}

      rot.x = 2 * S.s[j].x;
    } else {
      std::uniform_int_distribution<int> dist(0, 1);
      int j = dist(r);
      if (j) {
        rot.x = 2 * S.s[i0].x + 1;
      } else {
        rot.x = 2 * S.s[i0].x - 1;
      }
    }

    odd_run = !odd_run;

    return rot;
  };

  // initailze the measurement object
  simple_measurement A(S);

  // initialize the random number generator
  auto seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
  std::mt19937 rng{seed};

  // run wolff N times
  S.run_wolff(N, gen_R_IH, A, rng);

  // print the result of our measurements
  std::cout << "Wolff complete!\nThe average energy per site was " << A.avgE() / S.nv
    << ".\nThe average magnetization per site was " << A.avgM() / S.nv
    << ".\nThe average cluster size per site was " << A.avgC() / S.nv << ".\n";

  // exit
  return 0;
}

