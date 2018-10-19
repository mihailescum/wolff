
#include <getopt.h>
#include <iostream>
#include <chrono>

#include <GL/glut.h>

#include <wolff/models/ising.hpp>
#include <wolff/finite_states.hpp>

#include <wolff.hpp>

using namespace wolff;

class draw_ising : public measurement<ising_t, ising_t> {
  private:
    unsigned int frame_skip;
    v_t C;
  public:
    draw_ising(const system<ising_t, ising_t>& S, unsigned int window_size, unsigned int frame_skip, int argc, char *argv[]) : frame_skip(frame_skip){
      glutInit(&argc, argv);
      glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);
      glutInitWindowSize(window_size, window_size);
      glutCreateWindow("wolff");
      glClearColor(0.0,0.0,0.0,0.0);
      glMatrixMode(GL_PROJECTION);
      glLoadIdentity();
      gluOrtho2D(0.0, S.G.L, 0.0, S.G.L);
    }

    void pre_cluster(N_t, N_t, const system<ising_t, ising_t>& S, v_t, const ising_t&) {
      glClear(GL_COLOR_BUFFER_BIT);
      for (v_t i = 0; i < pow(S.G.L, 2); i++) {
        if (S.s[i].x == S.s0.x) {
          glColor3f(0.0, 0.0, 0.0);
        } else {
          glColor3f(1.0, 1.0, 1.0);
        }
        glRecti(i / S.G.L, i % S.G.L, (i / S.G.L) + 1, (i % S.G.L) + 1);
      }
      glFlush();
      C = 0;
    }

    void plain_bond_visited(const system<ising_t, ising_t>&, v_t, const ising_t&, v_t, double dE) {}

    void ghost_bond_visited(const system<ising_t, ising_t>&, v_t, const ising_t& s_old, const ising_t& s_new, double dE) {}

    void plain_site_transformed(const system<ising_t, ising_t>& S, v_t i, const ising_t&) {
      glColor3f(1.0, 0.0, 0.0);
      glRecti(i / S.G.L, i % S.G.L, (i / S.G.L) + 1, (i % S.G.L) + 1);
      C++;
      if (C % frame_skip == 0) {
        glFlush();
      }
    }

    void ghost_site_transformed(const system<ising_t, ising_t>&, const ising_t&) {}

    void post_cluster(N_t, N_t, const system<ising_t, ising_t>&) {}
};

int main(int argc, char *argv[]) {

  // set defaults
  N_t N = (N_t)1e4;
  D_t D = 2;
  L_t L = 128;
  double T = 2.26918531421;
  double H = 0.0;
  unsigned int window_size = 512;
  unsigned int frame_skip = 1;

  int opt;

  // take command line arguments
  while ((opt = getopt(argc, argv, "N:D:L:T:H:w:f:")) != -1) {
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
    case 'w':
      window_size = atoi(optarg);
      break;
    case 'f':
      frame_skip = atoi(optarg);
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

  // define function that generates self-inverse rotations
  std::function <ising_t(std::mt19937&, const system<ising_t, ising_t>&, v_t)> gen_R = [] (std::mt19937&, const system<ising_t, ising_t>&, v_t) -> ising_t {
    return ising_t(true);
  };

  // initailze the measurement object
  draw_ising A(S, window_size, frame_skip, argc, argv);

  // initialize the random number generator
  auto seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
  std::mt19937 rng{seed};

  // run wolff N times
  S.run_wolff(N, gen_R, A, rng);

  // exit
  return 0;
}

