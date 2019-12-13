
#include <getopt.h>
#include <iostream>
#include <chrono>

#include <wolff_models/height.hpp>
#include <wolff_models/dihedral_inf.hpp>
#include <wolff.hpp>

#include "simple_measurement.hpp"

int main(int argc, char *argv[]) {

  // set defaults
  unsigned N = (unsigned)1e4;
  unsigned D = 2;
  unsigned L = 128;
  double T = 0.8;
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
  std::function <double(const height_t<double>&, const height_t<double>&)> Z = [] (const height_t<double>& s1, const height_t<double>& s2) -> double {
    return - pow(s1.x - s2.x, 2);
  };

  // define the spin-field coupling
  std::function <double(const height_t<double>&)> B = [&] (const height_t<double>& s) -> double {
    return - H * pow(s.x, 2);
  };

  // initialize the lattice
  graph<> G(D, L);

  // initialize the system
  wolff::system<dihedral_inf_t<double>, height_t<double>, graph<>> S(G, T, Z, B);

  bool odd_run = false;

  std::function <dihedral_inf_t<double>(std::mt19937&, const wolff::system<dihedral_inf_t<double>, height_t<double>, graph<>>&, const graph<>::vertex&)> gen_R_IH = [&](std::mt19937& r, const wolff::system<dihedral_inf_t<double>, height_t<double>, graph<>>& S, const graph<>::vertex& v) -> dihedral_inf_t<double> {
    dihedral_inf_t<double> rot;
    rot.is_reflection = true;

    if (odd_run) {
      std::uniform_int_distribution<unsigned> dist(0, S.nv - 2);
      unsigned j = v.ind;

      //while (S.s[j].x == S.s[i0].x) {
        unsigned tmp = dist(r);

        if (tmp < v.ind) {
          j = tmp;
        } else {
          j = tmp + 1;
        }
      //}

      rot.x = 2 * S.s[j].x;
    } else {
      std::normal_distribution<double> dist(0.0,1.0);
      rot.x = 2 * S.s[v.ind].x + dist(r);
    }

    odd_run = !odd_run;

    return rot;
  };

  // initailze the measurement object
  simple_measurement A(S);

  // initialize the random number generator
  auto seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
  std::mt19937 rng(seed);

  // run wolff N times
  S.run_wolff(N, gen_R_IH, A, rng);

  // print the result of our measurements
  std::cout << "Wolff complete!\nThe average energy per site was " << A.avgE() / S.nv
    << ".\nThe average magnetization per site was " << A.avgM() / S.nv
    << ".\nThe average cluster size per site was " << A.avgC() / S.nv << ".\n";

  // exit
  return 0;
}

