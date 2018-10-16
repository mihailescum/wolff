
#include <getopt.h>
#include <iostream>

#include "include/randutils/randutils.hpp"

#include <wolff.hpp>

// define your R_t and X_t. here both are the same, ising_t
class ising_t {
  public:
    bool x;

    // both X_t and R_t need a default constructor (and destructor, if relevant)
    ising_t() : x(false) {}
    ising_t(bool x) : x(x) {}

    // R_t needs the member functions
    //   X_t act(const X_t& s) const {}
    //   R_t act(const R_t& s) const {}
    // to define action on both spins and other transformations
    ising_t act(const ising_t& s) const {
      if (x) {
        return ising_t(!(s.x));
      } else {
        return ising_t(s.x);
      }
    }

    // R_t needs the member functions
    //   X_t act_inverse(const X_t& s) const {}
    //   R_t act_inverse(const R_t& s) const {}
    // to define action of its inverse on both spins and other transformations
    ising_t act_inverse(const ising_t& s) const {
      return this->act(s);
    }
};

// define how measurements should be taken by importing the interface wolff_measurement<R_t, X_t>
class ising_measurements : public wolff_measurement<ising_t, ising_t> {
  private:
    count_t n;

    double E;
    int M;
    v_t S;

    double totalE;
    double totalM;
    double totalS;

  public:
    ising_measurements(D_t D, L_t L, double H) {
      n = 0;
      M = (int)pow(L, D);
      E = -D * pow(L, D) - H * pow(L, D);

      totalE = 0;
      totalM = 0;
      totalS = 0;
    }

    void pre_cluster(const state_t<ising_t, ising_t>& s, count_t step, count_t N, v_t v0, const ising_t& R) {
      S = 0;
    }

    void plain_bond_added(v_t v, const ising_t& s_old, const ising_t& s_new, v_t vn, const ising_t& sn, double dE) {
      E += dE;
    }

    void ghost_bond_added(v_t v, const ising_t& s_old, const ising_t& s_new, double dE) {
      E += dE;
      
      if (s_old.x) {
        M++;
      } else {
        M--;
      }

      if (s_new.x) {
        M--;
      } else {
        M++;
      }
    }

    void plain_site_transformed(v_t v, const ising_t& s_old, const ising_t& s_new) {
      S++;
    }

    void ghost_site_transformed(const ising_t& R_old, const ising_t& R_new) {
    }

    void post_cluster(const state_t<ising_t, ising_t>& s, count_t step, count_t N) {
      totalE += E;
      totalM += M;
      totalS += S;
      n++;
    }

    double avgE() {
      return totalE / n;
    }

    double avgM() {
      return totalM / n;
    }

    double avgS() {
      return totalS / n;
    }
};

int main(int argc, char *argv[]) {

  // set defaults
  count_t N = (count_t)1e4;
  D_t D = 2;
  L_t L = 128;
  double T = 2.26918531421;
  double H = 0.0;

  int opt;

  // take command line arguments
  while ((opt = getopt(argc, argv, "N:D:L:T:H:")) != -1) {
    switch (opt) {
    case 'N': // number of steps
      N = (count_t)atof(optarg);
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
    if (s1.x == s2.x) {
      return 1.0;
    } else {
      return -1.0;
    }
  };

  // define the spin-field coupling
  std::function <double(const ising_t&)> B = [=] (const ising_t& s) -> double {
    if (s.x) {
      return -H;
    } else {
      return H;
    }
  };

  // initialize the system
  state_t<ising_t, ising_t> s(D, L, T, Z, B);

  // initialize the random number generator
  randutils::auto_seed_128 seeds;
  std::mt19937 rng{seeds};

  // define function that generates self-inverse rotations
  std::function <ising_t(std::mt19937&, const ising_t&)> gen_R = [] (std::mt19937&, const ising_t& s) -> ising_t {
    return ising_t(true);
  };

  // initailze the measurement object
  ising_measurements m(D, L, H);

  // run wolff N times
  wolff<ising_t, ising_t>(N, s, gen_R, m, rng);

  // print the result of our measurements
  std::cout << "Wolff complete!\nThe average energy per site was " << m.avgE() / s.nv
    << ".\nThe average magnetization per site was " << m.avgM() / s.nv
    << ".\nThe average cluster size per site was " << m.avgS() / s.nv << ".\n";

  // exit
  return 0;
}

