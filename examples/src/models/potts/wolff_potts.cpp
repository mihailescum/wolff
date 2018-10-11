
#include <getopt.h>
#include <stdio.h>

#ifdef HAVE_GLUT
#include <GL/glut.h>
#endif

// include your group and spin space
#include "symmetric.hpp"
#include "potts.hpp"

// hack to speed things up considerably
#define N_STATES POTTSQ
#include <wolff/finite_states.hpp>

// include wolff.h
#include <measure.hpp>
#include <colors.h>
#include <randutils/randutils.hpp>
#include <wolff.hpp>

typedef state_t <symmetric_t<POTTSQ>, potts_t<POTTSQ>> sim_t;

int main(int argc, char *argv[]) {

  count_t N = (count_t)1e4;

  D_t D = 2;
  L_t L = 128;
  double T = 2.26918531421;
  double *H_vec = (double *)calloc(MAX_Q, sizeof(double));

  bool silent = false;
  bool draw = false;
  bool N_is_sweeps = false;
  unsigned int window_size = 512;

  // don't measure anything by default
  unsigned char measurement_flags = 0;

  int opt;
  q_t H_ind = 0;

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
    case 'H': // external field. nth call couples to state n
      H_vec[H_ind] = atof(optarg);
      H_ind++;
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
  std::function <double(const potts_t<POTTSQ>&, const potts_t<POTTSQ>&)> Z = [] (const potts_t<POTTSQ>& s1, const potts_t<POTTSQ>& s2) -> double {
    if (s1.x == s2.x) {
      return 1.0;
    } else {
      return 0.0;
    }
  };

  // define spin-field coupling
  std::function <double(const potts_t<POTTSQ> &)> B = [=] (const potts_t<POTTSQ>& s) -> double {
    return H_vec[s.x];
  };

  // initialize state object
  state_t <symmetric_t<POTTSQ>, potts_t<POTTSQ>> s(D, L, T, Z, B);

  // define function that generates self-inverse rotations
  std::function <symmetric_t<POTTSQ>(std::mt19937&, potts_t<POTTSQ>)> gen_R = [] (std::mt19937& r, potts_t<POTTSQ> v) -> symmetric_t<POTTSQ> {
    symmetric_t<POTTSQ> rot;

    std::uniform_int_distribution<q_t> dist(0, POTTSQ - 1);
    q_t j = dist(r);
    q_t swap_v;
    if (j < v.x) {
      swap_v = j;
    } else {
      swap_v = j + 1;
    }

    rot[v.x] = swap_v;
    rot[swap_v] = v.x;

    return rot;
  };

  FILE **outfiles = measure_setup_files(measurement_flags, timestamp);

  std::function <void(const sim_t&)> other_f;
  uint64_t sum_of_clusterSize = 0;

  if (N_is_sweeps) {
    other_f = [&] (const sim_t& s) {
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

    other_f = [] (const sim_t& s) {
      glClear(GL_COLOR_BUFFER_BIT);
      for (v_t i = 0; i < pow(s.L, 2); i++) {
        potts_t<POTTSQ> tmp_s = s.R.act_inverse(s.spins[i]);
        glColor3f(hue_to_R(tmp_s.x * 2 * M_PI / POTTSQ), hue_to_G(tmp_s.x * 2 * M_PI / POTTSQ), hue_to_B(tmp_s.x * 2 * M_PI / POTTSQ));
        glRecti(i / s.L, i % s.L, (i / s.L) + 1, (i % s.L) + 1);
      }
      glFlush();
    };
#endif
  } else {
    other_f = [] (const sim_t& s) {};
  }

  std::function <void(const sim_t&)> measurements = measure_function_write_files(measurement_flags, outfiles, other_f);

  // add line to metadata file with run info
  {
    FILE *outfile_info = fopen("wolff_metadata.txt", "a");

    fprintf(outfile_info, "<| \"ID\" -> %lu, \"MODEL\" -> \"POTTS\", \"q\" -> %d, \"D\" -> %" PRID ", \"L\" -> %" PRIL ", \"NV\" -> %" PRIv ", \"NE\" -> %" PRIv ", \"T\" -> %.15f, \"H\" -> {", timestamp, POTTSQ, s.D, s.L, s.nv, s.ne, T);

    for (q_t i = 0; i < POTTSQ; i++) {
      fprintf(outfile_info, "%.15f", H_vec[i]);
      if (i < POTTSQ - 1) {
        fprintf(outfile_info, ", ");
      }
    }

    fprintf(outfile_info, "} |>\n");

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

  // free the random number generator
  free(H_vec);
  measure_free_files(measurement_flags, outfiles);

  return 0;

}

