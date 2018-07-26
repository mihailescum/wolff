
#include <getopt.h>

#ifdef HAVE_GLUT
#include <GL/glut.h>
#endif

// include your group and spin space
#include <dihedral_inf.h>
#include <height.h>

// include wolff.h
#include <rand.h>
#include <wolff.h>

typedef state_t <dihedral_inf_t<double>, height_t<double>> sim_t;

int main(int argc, char *argv[]) {

  count_t N = (count_t)1e4;

  D_t D = 2;
  L_t L = 128;
  double T = 2.26918531421;
  double H = 0;

  bool silent = false;
  bool draw = false;
  unsigned int window_size = 512;
  double epsilon = 1;

  int opt;

  while ((opt = getopt(argc, argv, "N:D:L:T:H:sdw:e:")) != -1) {
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
      H = atof(optarg);
      break;
    case 'e': // external field. nth call couples to state n
      epsilon = atof(optarg);
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
  std::function <double(const height_t<double>&, const height_t<double>&)> Z = [] (const height_t<double>& h1, const height_t<double>& h2) -> double {
    return -pow(h1.x - h2.x, 2);
  };

  // define spin-field coupling
  std::function <double(height_t<double>)> B = [=] (height_t<double> h) -> double {
    return -H * pow(h.x, 2);;
  };

  // initialize state object
  sim_t s(D, L, T, Z, B);

  // define function that generates self-inverse rotations
  std::function <dihedral_inf_t<double>(gsl_rng *, height_t<double>)> gen_R = [=] (gsl_rng *r, height_t<double> h) -> dihedral_inf_t<double> {
    dihedral_inf_t<double> rot;
    rot.is_reflection = true;

    double amount = epsilon * gsl_ran_ugaussian(r);

    rot.x = 2 * h.x + amount;

    return rot;
  };

  // define function that updates any number of measurements
  std::function <void(const sim_t&)> measurement;

  double average_M = 0;
  if (!draw) {
    // a very simple example: measure the average magnetization
    measurement = [&] (const sim_t& s) {
      average_M += (double)s.M / (double)N / (double)s.nv;
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

    measurement = [&] (const sim_t& s) {
      average_M += (double)s.M / (double)N / (double)s.nv;
      glClear(GL_COLOR_BUFFER_BIT);
      double max_h = INT64_MIN;
      double min_h = INT64_MAX;
      for (v_t i = 0; i < pow(L, 2); i++) {
        double cur_h = (s.R.act_inverse(s.spins[i])).x;
        if (cur_h < min_h) {
          min_h = cur_h;
        }
        if (cur_h > max_h) {
          max_h = cur_h;
        }
      }

      for (v_t i = 0; i < pow(L, 2); i++) {
        double cur_h = (s.R.act_inverse(s.spins[i])).x;
        double mag = ((double)(cur_h - min_h)) / ((double)(max_h - min_h));
        glColor3f(mag, mag, mag);
        glRecti(i / L, i % L, (i / L) + 1, (i % L) + 1);
      }
      glFlush();
    };
#endif
  }

  // run wolff for N cluster flips
  wolff(N, s, gen_R, measurement, r, silent);

  // tell us what we found!
  printf("%" PRIcount " DGM runs completed. D = %" PRID ", L = %" PRIL ", T = %g, H = %g, <M> = %g\n", N, D, L, T, H, average_M);

  // free the random number generator
  gsl_rng_free(r);

  if (draw) {
  }

  return 0;

}

