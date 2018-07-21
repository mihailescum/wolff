
#include <getopt.h>
#ifdef HAVE_GLUT
#include <GL/glut.h>
#endif

#include <wolff.h>
#include <correlation.h>
#include <measure.h>

typedef orthogonal_t <N_COMP, double> orthogonal_R_t;
typedef vector_t <N_COMP, double> vector_R_t;
typedef state_t <orthogonal_R_t, vector_R_t> On_t;

// angle from the x-axis of a two-vector
double theta(vector_R_t v) {
  double x = v.x[0];
  double y = v.x[1];

  double val = atan(y / x);

  if (x < 0.0 && y > 0.0) {
    return M_PI + val;
  } else if ( x < 0.0 && y < 0.0 ) {
    return - M_PI + val;
  } else {
    return val;
  }
}

double H_modulated(vector_R_t v, int order, double mag) {
  return mag * cos(order * theta(v));
}

double hue_to_R(double theta) {
  if (((M_PI / 3 <= theta) && (theta < 2 * M_PI / 3)) || ((4 * M_PI / 3 <= theta) && (theta < 5 * M_PI / 3))) {
    return 1.0 - fabs(fmod(theta / (2 * M_PI / 6), 2) - 1.0);
  } else if (((0 <= theta) && (theta < M_PI / 3)) || ((5 * M_PI / 3 <= theta) && (theta <= 2 * M_PI))) {
    return 1.0;
  } else {
    return 0.0;
  }
}

double hue_to_G(double theta) {
  if (((0 <= theta) && (theta < M_PI / 3)) || ((M_PI <= theta) && (theta < 4 * M_PI / 3))) {
    return 1.0 - fabs(fmod(theta / (2 * M_PI / 6), 2) - 1.0);
  } else if (((M_PI / 3 <= theta) && (theta < 2 * M_PI / 3)) || ((2 * M_PI / 3 <= theta) && (theta < M_PI))) {
    return 1.0;
  } else {
    return 0.0;
  }
}

double hue_to_B(double theta) {
  if (((2 * M_PI / 3 <= theta) && (theta < M_PI)) || ((5 * M_PI / 3 <= theta) && (theta <= 2 * M_PI))) {
    return 1.0 - fabs(fmod(theta / (2 * M_PI / 6), 2) - 1.0);
  } else if (((M_PI <= theta) && (theta < 4 * M_PI / 3)) || ((4 * M_PI / 3 <= theta) && (theta < 5 * M_PI / 3))) {
    return 1.0;
  } else {
    return 0.0;
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
  int order = 2;

  int opt;
  q_t J_ind = 0;
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

  std::function <orthogonal_R_t(gsl_rng *, const On_t *)> gen_R;

  if (use_pert) {
    gen_R = std::bind(generate_rotation_perturbation <N_COMP>, std::placeholders::_1, std::placeholders::_2, epsilon);
    pert_type = "PERTURB";
  } else {
    gen_R = generate_rotation_uniform <N_COMP>;
    pert_type = "UNIFORM";
  }

  FILE *outfile_info = fopen("wolff_metadata.txt", "a");

  fprintf(outfile_info, "<| \"ID\" -> %lu, \"MODEL\" -> \"%s\", \"q\" -> %d, \"D\" -> %" PRID ", \"L\" -> %" PRIL ", \"NV\" -> %" PRIv ", \"NE\" -> %" PRIv ", \"T\" -> %.15f, \"FIELD_TYPE\" -> ", timestamp, ON_strings[N_COMP], N_COMP, D, L, (v_t)pow(L, D), D * (v_t)pow(L, D), T);
  if (modulated_field) {
    fprintf(outfile_info, "\"MODULATED\", \"ORDER\" -> %d, \"H\" -> %.15f, ", order, H_vec[0]);
  } else {
    fprintf(outfile_info, "\"VECTOR\", \"H\" -> {");
    for (q_t i = 0; i < N_COMP; i++) {
      fprintf(outfile_info, "%.15f", H_vec[i]);
      if (i < N_COMP - 1) {
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
        vector_R_t v_tmp = act_inverse(s->R, s->spins[i]);
        double thetai = fmod(2 * M_PI + theta(v_tmp), 2 * M_PI);
        free_spin(v_tmp);
        double saturation = 0.7;
        double value = 0.9;
        double chroma = saturation * value;
        glColor3f(chroma * hue_to_R(thetai) + (value - chroma), chroma * hue_to_G(thetai) + (value - chroma), chroma * hue_to_B(thetai) + (value - chroma));
        glRecti(i / L, i % L, (i / L) + 1, (i % L) + 1);
      }
      glFlush();
    };
#endif
  } else {
    other_f = [] (const On_t *s) {};
  }

  std::function <void(const On_t *)> measurements = measure_function_write_files(measurement_flags, outfiles, other_f);

  std::function <double(vector_R_t)> H;

  if (modulated_field) {
    H = std::bind(H_modulated, std::placeholders::_1, order, H_vec[0]);
  } else {
    H = std::bind(H_vector <N_COMP, double>, std::placeholders::_1, H_vec);
  }

  // initialize random number generator
  gsl_rng *r = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(r, rand_seed());

  state_t <orthogonal_R_t, vector_R_t> s(D, L, T, dot <N_COMP, double>, H);

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

