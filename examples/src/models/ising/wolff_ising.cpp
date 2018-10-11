
#include <getopt.h>
#include <stdio.h>

// if you have GLUT installed, you can see graphics!
#ifdef HAVE_GLUT
#include <GL/glut.h>
#endif

// include your group and spin space
#include "z2.hpp"
#include "ising.hpp"

// finite_states.h can be included for spin types that have special variables
// defined, and it causes wolff execution to use precomputed bond probabilities
#include <wolff/finite_states.hpp>

#include <randutils/randutils.hpp>

// measure.hpp contains useful functions for saving timeseries to files
#include <measure.hpp>

// include wolff.hpp
#include <wolff.hpp>

int main(int argc, char *argv[]) {

  count_t N = (count_t)1e4;

  D_t D = 2;
  L_t L = 128;
  double T = 2.26918531421;
  double H = 0.0;

  bool silent = false;
  bool draw = false;
  bool N_is_sweeps = false;
  unsigned int window_size = 512;

  // don't measure anything by default
  unsigned char measurement_flags = 0;

  int opt;

  while ((opt = getopt(argc, argv, "N:D:L:T:H:sdw:M:S")) != -1) {
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
    case 's': // don't print anything during simulation. speeds up slightly
      silent = true;
      break;
    case 'S':
      N_is_sweeps = true;
      break;
    case 'd':
#ifdef HAVE_GLUT
      draw = true;
      break;
#else
      printf("You didn't compile this with the glut library installed!\n");
      exit(EXIT_FAILURE);
#endif
    case 'w':
      window_size = atoi(optarg);
      break;
    case 'M':
      measurement_flags ^= 1 << atoi(optarg);
      break;
    default:
      exit(EXIT_FAILURE);
    }
  }

  // get nanosecond timestamp for unique run id
  unsigned long timestamp;

  {
    struct timespec spec;
    clock_gettime(CLOCK_REALTIME, &spec);
    timestamp = spec.tv_sec*1000000000LL + spec.tv_nsec;
  }

  // initialize random number generator
  randutils::auto_seed_128 seeds;
  std::mt19937 rng{seeds};

  // define spin-spin coupling
  std::function <double(const ising_t&, const ising_t&)> Z = [] (const ising_t& s1, const ising_t& s2) -> double {
    if (s1.x == s2.x) {
      return 1.0;
    } else {
      return -1.0;
    }
  };

  // define spin-field coupling
  std::function <double(const ising_t&)> B = [=] (const ising_t& s) -> double {
    if (s.x) {
      return -H;
    } else {
      return H;
    }
  };

  // initialize state object
#ifndef NOFIELD
  state_t <z2_t, ising_t> s(D, L, T, Z, B);
#else
  state_t <z2_t, ising_t> s(D, L, T, Z);
#endif

  // define function that generates self-inverse rotations
  std::function <z2_t(std::mt19937&, ising_t)> gen_R = [] (std::mt19937&, const ising_t& s) -> z2_t {
    return z2_t(true);
  };

  FILE **outfiles = measure_setup_files(measurement_flags, timestamp);

  std::function <void(const state_t<z2_t, ising_t>&)> other_f;
  uint64_t sum_of_clusterSize = 0;

  if (N_is_sweeps) {
    other_f = [&] (const state_t<z2_t, ising_t>& s) {
      sum_of_clusterSize += s.last_cluster_size;
    };
  } else if (draw) {
#ifdef HAVE_GLUT
    // initialize glut
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);
    glutInitWindowSize(window_size, window_size);
    glutCreateWindow("wolff");
    glClearColor(0.0,0.0,0.0,0.0);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluOrtho2D(0.0, L, 0.0, L);

    other_f = [] (const state_t <z2_t, ising_t>& s) {
      glClear(GL_COLOR_BUFFER_BIT);
      for (v_t i = 0; i < pow(s.L, 2); i++) {
        if (s.spins[i].x == s.R.x) {
          glColor3f(0.0, 0.0, 0.0);
        } else {
          glColor3f(1.0, 1.0, 1.0);
        }
        glRecti(i / s.L, i % s.L, (i / s.L) + 1, (i % s.L) + 1);
      }
      glFlush();
    };
#endif
  } else {
    other_f = [] (const state_t<z2_t, ising_t>& s) {};
  }

  std::function <void(const state_t<z2_t, ising_t>&)> measurements = measure_function_write_files(measurement_flags, outfiles, other_f);

  // add line to metadata file with run info
  {
    FILE *outfile_info = fopen("wolff_metadata.txt", "a");

    fprintf(outfile_info, "<| \"ID\" -> %lu, \"MODEL\" -> \"ISING\", \"q\" -> 2, \"D\" -> %" PRID ", \"L\" -> %" PRIL ", \"NV\" -> %" PRIv ", \"NE\" -> %" PRIv ", \"T\" -> %.15f, \"H\" -> %.15f |>\n", timestamp, s.D, s.L, s.nv, s.ne, T, H);

    fclose(outfile_info);
  }

  // run wolff for N cluster flips
  if (N_is_sweeps) {
    count_t N_rounds = 0;
    printf("\n");
    while (sum_of_clusterSize < N * s.nv) {
      printf("\033[F\033[J\033[F\033[JWOLFF: sweep %" PRIu64 " / %" PRIu64 ": E = %.2f, S = %" PRIv "\n", (count_t)((double)sum_of_clusterSize / (double)s.nv), N, s.E, s.last_cluster_size);
      wolff(N, s, gen_R, measurements, rng, silent);
      N_rounds++;
    }
    printf("\033[F\033[J\033[F\033[JWOLFF: sweep %" PRIu64 " / %" PRIu64 ": E = %.2f, S = %" PRIv "\n\n", (count_t)((double)sum_of_clusterSize / (double)s.nv), N, s.E, s.last_cluster_size);
  } else {
    wolff(N, s, gen_R, measurements, rng, silent);
  }

  measure_free_files(measurement_flags, outfiles);

  return 0;

}

