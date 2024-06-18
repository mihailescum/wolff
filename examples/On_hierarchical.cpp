
#include <getopt.h>
#include <iostream>
#include <chrono>

#include <wolff.hpp>
#include <wolff_models/vector.hpp>
#include <wolff_models/orthogonal.hpp>
#include <wolff_graphs/hierarchical.hpp>

#include "export_measurement.hpp"

typedef wolff::hierarchical_graph<> graph_type;

int main(int argc, char *argv[])
{

  // set defaults
  unsigned N = (unsigned)1e2;
  unsigned D = WOLFF_D;
  unsigned R = 1;
  unsigned L = 128;
  double T = 0.8;
  wolff::vector_t<WOLFF_N, double> H;
  H.fill(0.0);
  unsigned Hi = 0;

  int opt;

  // take command line arguments
  while ((opt = getopt(argc, argv, "N:D:L:T:H:R:")) != -1)
  {
    switch (opt)
    {
    case 'N': // number of steps
      N = (unsigned)atof(optarg);
      break;
    case 'L': // linear size
      L = atoi(optarg);
      break;
    case 'R': // linear size
      R = atoi(optarg);
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
  std::function<double(const graph_type::halfedge &, const wolff::vector_t<WOLFF_N, double> &, const wolff::vector_t<WOLFF_N, double> &)> Z = [](const graph_type::halfedge &edge, const wolff::vector_t<WOLFF_N, double> &s1, const wolff::vector_t<WOLFF_N, double> &s2) -> double
  {
    return edge.prop * s1 * s2;
  };

  // define the spin-field coupling
  std::function<double(const wolff::vector_t<WOLFF_N, double> &)> B = [&](const wolff::vector_t<WOLFF_N, double> &s) -> double
  {
    return H * s;
  };

  // initialize the lattice
  graph_type G(D, L, R);
  G.init();

  // initialize the system
  wolff::system<wolff::orthogonal_t<WOLFF_N, double>, wolff::vector_t<WOLFF_N, double>, graph_type> S(G, T, Z, B);

  std::function<wolff::orthogonal_t<WOLFF_N, double>(std::mt19937 &, const wolff::system<wolff::orthogonal_t<WOLFF_N, double>, wolff::vector_t<WOLFF_N, double>, graph_type> &, const graph_type::vertex)> gen_R = wolff::generate_rotation_uniform<WOLFF_N, graph_type>;

  // initailze the measurement object
  export_measurement<wolff::orthogonal_t<WOLFF_N, double>, wolff::vector_t<WOLFF_N, double>, graph_type> M(S);

  // initialize the random number generator
  auto seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
  std::mt19937 rng(seed);

  // run wolff N times
  S.run_wolff(N, gen_R, M, rng);

  M.export_results("hierarchical_out.txt");

  // exit
  return 0;
}
