
#include <getopt.h>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <chrono>

#include <wolff_models/vector.hpp>
#include <wolff_models/orthogonal.hpp>

#include "simple_measurement.hpp"

int main(int argc, char *argv[])
{

  // set defaults
  unsigned N = (unsigned)1e2;
  unsigned D = 2;
  unsigned L = 128;
  double T = 0.8;
  vector_t<WOLFF_N, double> H;
  H.fill(0.0);
  unsigned Hi = 0;

  int opt;

  // take command line arguments
  while ((opt = getopt(argc, argv, "N:D:L:T:H:")) != -1)
  {
    switch (opt)
    {
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
  std::function<double(const vector_t<WOLFF_N, double> &, const vector_t<WOLFF_N, double> &)> Z = [](const vector_t<WOLFF_N, double> &s1, const vector_t<WOLFF_N, double> &s2) -> double
  {
    return s1 * s2;
  };

  // define the spin-field coupling
  std::function<double(const vector_t<WOLFF_N, double> &)> B = [&](const vector_t<WOLFF_N, double> &s) -> double
  {
    return H * s;
  };

  // initialize the lattice
  graph<> G(D, L);

  // initialize the system
  wolff::system<orthogonal_t<WOLFF_N, double>, vector_t<WOLFF_N, double>> S(G, T, Z, B);

  std::function<orthogonal_t<WOLFF_N, double>(std::mt19937 &, const wolff::system<orthogonal_t<WOLFF_N, double>, vector_t<WOLFF_N, double>, graph<>> &, const graph<>::vertex)> gen_R = generate_rotation_uniform<WOLFF_N, graph<>>;

  // initailze the measurement object
  simple_measurement A(S);

  // initialize the random number generator
  auto seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
  std::mt19937 rng(seed);

  // run wolff N times
  S.run_wolff(N, gen_R, A, rng);

  // print the result of our measurements
  std::cout << "Wolff complete!\nThe average energy per site was " << A.avgE() / S.nv
            << ".\nThe average magnetization per site was " << A.avgM() / S.nv
            << ".\nThe average cluster size per site was " << A.avgC() / S.nv << ".\n";

  std::stringstream result_array;
  for (const auto &s : S.s)
  {
    for (const auto &e : s)
    {
      result_array << std::setprecision(4) << e << " ";
    }
    // result_array << "\n";
  }

  std::ofstream outFile("on_result.txt");
  outFile << result_array.rdbuf();
  outFile.close();

  // exit
  return 0;
}
