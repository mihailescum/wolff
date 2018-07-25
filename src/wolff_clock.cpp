
#include <getopt.h>

#ifdef HAVE_GLUT
#include <GL/glut.h>
#endif

// include your group and spin space
#include <dihedral.h>
#include <potts.h>
#include <colors.h>

// hack to speed things up considerably
#define N_STATES POTTSQ
#include <finite_states.h>

// include wolff.h
#include <rand.h>
#include <wolff.h>

typedef state_t <dihedral_t<q_t,POTTSQ>, potts_t<POTTSQ>> sim_t;

int main(int argc, char *argv[]) {

  count_t N = (count_t)1e4;

  D_t D = 2;
  L_t L = 128;
  double T = 2.26918531421;
  double *H_vec = (double *)calloc(MAX_Q, sizeof(double));

  bool silent = false;
  bool draw = false;
  unsigned int window_size = 512;

  int opt;
  q_t H_ind = 0;

  while ((opt = getopt(argc, argv, "N:D:L:T:H:sdw:")) != -1) {
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
    default:
      exit(EXIT_FAILURE);
    }
  }

  // initialize random number generator
  gsl_rng *r = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(r, rand_seed());

  // define spin-spin coupling
  std::function <double(potts_t<POTTSQ>, potts_t<POTTSQ>)> Z = [] (potts_t<POTTSQ> s1, potts_t<POTTSQ> s2) -> double {
    return cos(2 * M_PI * (double)(s1.x + POTTSQ - s2.x) / (double)POTTSQ);
  };

  // define spin-field coupling
  std::function <double(potts_t<POTTSQ>)> B = [=] (potts_t<POTTSQ> s) -> double {
    return H_vec[s.x];
  };

  // initialize state object
  state_t <dihedral_t<q_t,POTTSQ>, potts_t<POTTSQ>> s(D, L, T, Z, B);

  // define function that generates self-inverse rotations
  std::function <dihedral_t<q_t,POTTSQ>(gsl_rng *, potts_t<POTTSQ>)> gen_R = [] (gsl_rng *r, potts_t<POTTSQ> v) -> dihedral_t<q_t,POTTSQ> {
    dihedral_t<q_t,POTTSQ> rot;
    rot.is_reflection = true;
    q_t x = gsl_rng_uniform_int(r, POTTSQ - 1);
    rot.x = (2 * v.x + x + 1) % POTTSQ;

    return rot;
  };

  // define function that updates any number of measurements
  std::function <void(const sim_t *)> measurement;

  double average_M = 0;
  if (!draw) {
    // a very simple example: measure the average magnetization
    measurement = [&] (const sim_t *s) {
      average_M += (double)s->M.x[0] / (double)N / (double)s->nv;
    };
  } else {
    // a more complex example: measure the average magnetization, and draw the spin configuration to the screen

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

    measurement = [&] (const sim_t *s) {
      average_M += (double)s->M.x[0] / (double)N / (double)s->nv;
      glClear(GL_COLOR_BUFFER_BIT);
      for (v_t i = 0; i < pow(L, 2); i++) {
        potts_t<POTTSQ> tmp_s = act_inverse(s->R, s->spins[i]);
        glColor3f(hue_to_R(tmp_s.x * 2 * M_PI / POTTSQ), hue_to_G(tmp_s.x * 2 * M_PI / POTTSQ), hue_to_B(tmp_s.x * 2 * M_PI / POTTSQ));
        glRecti(i / L, i % L, (i / L) + 1, (i % L) + 1);
      }
      glFlush();
    };
#endif
  }

  // run wolff for N cluster flips
  wolff(N, &s, gen_R, measurement, r, silent);

  // tell us what we found!
  printf("%" PRIcount " %d-Potts runs completed. D = %" PRID ", L = %" PRIL ", T = %g, H = %g, <M> = %g\n", N, POTTSQ, D, L, T, H_vec[0], average_M);

  // free the random number generator
  gsl_rng_free(r);

  if (draw) {
  }

  return 0;

}

