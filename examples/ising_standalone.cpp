#include <iostream>
#include <chrono>

#include <wolff.hpp>

using namespace wolff;

class ising_t {
  public:
    int s;

    ising_t() : s(1) {};
    ising_t(int i) : s(i) {};

    ising_t act(const ising_t& x) const {
      return ising_t(s * x.s);
    }

    ising_t act_inverse(const ising_t& x) const {
      return this->act(x);
    }
};

class measure_clusters : public measurement<ising_t, ising_t> {
  private:
    v_t C;

  public:
    double Ctotal;

    measure_clusters() { Ctotal = 0; }

    void pre_cluster(N_t, N_t, const system<ising_t, ising_t>&, v_t, const ising_t&) { C = 0; }

    void plain_site_transformed(const system<ising_t, ising_t>&, v_t, const ising_t&) { C++; }

    void post_cluster(N_t, N_t, const system<ising_t, ising_t>&) { Ctotal += C; }
};

int main(int argc, char *argv[]) {
  // set defaults
  N_t N = (N_t)1e3;
  D_t D = 2;
  L_t L = 128;
  double T = 2.26918531421;
  double H = 0.01;

  // define the spin-spin coupling
  std::function <double(const ising_t&, const ising_t&)> Z =
    [](const ising_t& s1, const ising_t& s2) -> double {
      return (double)(s1.s * s2.s);
    };

  // define the spin-field coupling
  std::function <double(const ising_t&)> B =
    [=](const ising_t& s) -> double {
      return H * s.s;
    };

  // initialize the lattice
  graph G(D, L);

  // initialize the system
  system<ising_t, ising_t> S(G, T, Z, B);

  // define function that generates self-inverse rotations
  std::function <ising_t(std::mt19937&, const system<ising_t, ising_t>&, v_t)> gen_R =
    [] (std::mt19937&, const system<ising_t, ising_t>&, v_t) -> ising_t {
      return ising_t(-1);
    };

  // initailze the measurement object
  measure_clusters A;

  // initialize the random number generator
  auto seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
  std::mt19937 rng{seed};

  // run wolff N times
  S.run_wolff(N, gen_R, A, rng);

  // print results
  std::cout << "The average cluster size per site was " << (A.Ctotal / N) / S.nv << ".\n";

  // exit
  return 0;
}
