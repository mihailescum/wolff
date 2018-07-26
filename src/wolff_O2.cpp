
#include <getopt.h>
#include <stdio.h>

#ifdef HAVE_GLUT
#include <GL/glut.h>
#endif

#include <circle_group.h>
#include <angle.h>

#include <wolff.h>
#include <measure.h>
#include <colors.h>
#include <rand.h>

typedef circle_group_t orthogonal_R_t;
typedef angle_t vector_R_t;
typedef state_t <orthogonal_R_t, vector_R_t> On_t;

double H_modulated(vector_R_t v, int order, double mag) {
  return mag * cos(order * v.x);
}

double theta(double *v) {
  double x = v[0];
  double y = v[1];

  if (x == 0) { 
    if (y >= 0) {
      return M_PI / 2;
    } else {
      return - M_PI / 2;
    }
  } else {
    double val = atan(y / x);

    if (x < 0.0 && y > 0.0) {
      return M_PI + val;
    } else if ( x < 0.0 && y < 0.0 ) {
      return - M_PI + val;
    } else {
      return val;
    }
  }
}

int main(int argc, char *argv[]) {

  count_t N = (count_t)1e7;

  D_t D = 2;
  L_t L = 128;
  double T = 2.26918531421;
  double *H_vec = (double *)calloc(MAX_Q, sizeof(double));

  bool silent = false;
  bool use_pert = false;
  bool N_is_sweeps = false;
  bool draw = false;
  unsigned int window_size = 512;

  bool modulated_field = false;
  unsigned int order = 1;

  int opt;
  q_t H_ind = 0;
  double epsilon = 1;

//  unsigned char measurement_flags = measurement_energy | measurement_clusterSize;

  unsigned char measurement_flags = 0;

  while ((opt = getopt(argc, argv, "N:D:L:T:H:spe:mo:M:Sdw:")) != -1) {
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
    case 'p':
      use_pert = true;
      break;
    case 'e':
      epsilon = atof(optarg);
      break;
    case 'm':
      modulated_field = true;
      break;
    case 'M':
      measurement_flags ^= 1 << atoi(optarg);
      break;
    case 'o':
      order = atoi(optarg);
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
    default:
      exit(EXIT_FAILURE);
    }
  }

  unsigned long timestamp;

  {
    struct timespec spec;
    clock_gettime(CLOCK_REALTIME, &spec);
    timestamp = spec.tv_sec*1000000000LL + spec.tv_nsec;
  }

  const char *pert_type;

  std::function <orthogonal_R_t(gsl_rng *, vector_R_t)> gen_R;

  if (use_pert) {
    gen_R = [=] (gsl_rng *r, const angle_t& t) -> circle_group_t {
      circle_group_t rot;
      rot.is_reflection = true;

      unsigned int x = gsl_rng_uniform_int(r, order);
      double amount = epsilon * gsl_ran_ugaussian(r);

      rot.x = fmod(2 * M_PI * (1.0 + (double)x / (double)order + amount), 2 * M_PI);

      return rot;
    };
    pert_type = "PERTURB";
  } else {
    gen_R = [=] (gsl_rng *r, const angle_t& t) -> circle_group_t {
      circle_group_t rot;
      rot.is_reflection = true;
      rot.x = 2 * M_PI * gsl_rng_uniform(r);

      return rot;
    };
    pert_type = "UNIFORM";
  }

  FILE *outfile_info = fopen("wolff_metadata.txt", "a");

  fprintf(outfile_info, "<| \"ID\" -> %lu, \"MODEL\" -> \"%s\", \"q\" -> %d, \"D\" -> %" PRID ", \"L\" -> %" PRIL ", \"NV\" -> %" PRIv ", \"NE\" -> %" PRIv ", \"T\" -> %.15f, \"FIELD_TYPE\" -> ", timestamp, ON_strings[2], 2, D, L, (v_t)pow(L, D), D * (v_t)pow(L, D), T);
  if (modulated_field) {
    fprintf(outfile_info, "\"MODULATED\", \"ORDER\" -> %d, \"H\" -> %.15f, ", order, H_vec[0]);
  } else {
    fprintf(outfile_info, "\"VECTOR\", \"H\" -> {");
    for (q_t i = 0; i < 2; i++) {
      fprintf(outfile_info, "%.15f", H_vec[i]);
      if (i < 2 - 1) {
        fprintf(outfile_info, ", ");
      }
    }
    fprintf(outfile_info, "}, ");
  }

  fprintf(outfile_info, "\"GENERATOR\" -> \"%s\"", pert_type);

  if (use_pert) {
    fprintf(outfile_info, ", \"EPS\" -> %g", epsilon);
  }

  fprintf(outfile_info, " |>\n");

  fclose(outfile_info);

  FILE **outfiles = measure_setup_files(measurement_flags, timestamp);

  std::function <void(const On_t *)> other_f;
  uint64_t sum_of_clusterSize = 0;

  if (N_is_sweeps) {
    other_f = [&] (const On_t *s) {
      sum_of_clusterSize += s->last_cluster_size;
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

    other_f = [&] (const On_t *s) {
      glClear(GL_COLOR_BUFFER_BIT);
      for (v_t i = 0; i < pow(L, 2); i++) {
        vector_R_t v_tmp = s->R.act_inverse(s->spins[i]);
        double saturation = 0.7;
        double value = 0.9;
        double chroma = saturation * value;
        glColor3f(chroma * hue_to_R(v_tmp.x) + (value - chroma), chroma * hue_to_G(v_tmp.x) + (value - chroma), chroma * hue_to_B(v_tmp.x) + (value - chroma));
        glRecti(i / L, i % L, (i / L) + 1, (i % L) + 1);
      }
      glFlush();
    };
#endif
  } else {
    other_f = [] (const On_t *s) {};
  }

  std::function <void(const On_t *)> measurements = measure_function_write_files(measurement_flags, outfiles, other_f);

  std::function <double(const angle_t&, const angle_t&)> J = [] (const angle_t& t1, const angle_t& t2) -> double {
    return cos(t1.x - t2.x);
  };

  std::function <double(const angle_t &)> H;

  if (modulated_field) {
    H = [=] (const angle_t& t) -> double {
      return H_vec[0] * cos(order * t.x);
    };
  } else {
    double mag = 0;
    for (q_t i = 0; i < 2; i++) {
      mag += pow(H_vec[i], 2);
    }
    mag = sqrt(mag);
    double t0 = theta(H_vec);
    H = [=] (const angle_t& t) -> double {
      return mag * cos(t0 + t.x);
    };
  }

  // initialize random number generator
  gsl_rng *r = gsl_rng_alloc(gsl_rng_taus2);
  gsl_rng_set(r, rand_seed());

  state_t <orthogonal_R_t, vector_R_t> s(D, L, T, J, H);

  if (N_is_sweeps) {
    count_t N_rounds = 0;
    printf("\n");
    while (sum_of_clusterSize < N * s.nv) {
      printf("\033[F\033[J\033[F\033[JWOLFF: sweep %" PRIu64 " / %" PRIu64 ": E = %.2f, S = %" PRIv "\n", (count_t)((double)sum_of_clusterSize / (double)s.nv), N, s.E, s.last_cluster_size);
      wolff <orthogonal_R_t, vector_R_t> (N, &s, gen_R, measurements, r, silent);
      N_rounds++;
    }
    printf("\033[F\033[J\033[F\033[JWOLFF: sweep %" PRIu64 " / %" PRIu64 ": E = %.2f, S = %" PRIv "\n\n", (count_t)((double)sum_of_clusterSize / (double)s.nv), N, s.E, s.last_cluster_size);
  } else {
    wolff <orthogonal_R_t, vector_R_t> (N, &s, gen_R, measurements, r, silent);
  }

  measure_free_files(measurement_flags, outfiles);
  free(H_vec);
  gsl_rng_free(r);

  return 0;
}

