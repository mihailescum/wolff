
#include <getopt.h>

#include <wolff.h>
#include <correlation.h>

typedef state_t <orthogonal_t <2, double>, vector_t <2, double>> planar_t;

// angle from the x-axis of a two-vector
double theta(vector_t <2, double> v) {
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

double H_modulated(vector_t <2, double> v, int order, double mag) {
  return mag * cos(order * theta(v));
}

int main(int argc, char *argv[]) {

  count_t N = (count_t)1e7;

  D_t D = 2;
  L_t L = 128;
  double T = 2.26918531421;
  double *H_vec = (double *)calloc(MAX_Q, sizeof(double));

  bool silent = false;
  bool use_pert = false;

  bool modulated_field = false;
  int order = 2;

  int opt;
  q_t J_ind = 0;
  q_t H_ind = 0;
  double epsilon = 1;

  while ((opt = getopt(argc, argv, "N:q:D:L:T:J:H:spe:mo:")) != -1) {
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
    case 'o':
      order = atoi(optarg);
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

  std::function <orthogonal_t <2, double>(gsl_rng *, const planar_t *)> gen_R;

  if (use_pert) {
    gen_R = std::bind(generate_rotation_perturbation <2>, std::placeholders::_1, std::placeholders::_2, epsilon);
    pert_type = "PERTURB";
  } else {
    gen_R = generate_rotation_uniform <2>;
    pert_type = "UNIFORM";
  }


  FILE *outfile_info = fopen("wolff_metadata.txt", "a");

  fprintf(outfile_info, "<| \"ID\" -> %lu, \"MODEL\" -> \"PLANAR\", \"q\" -> 2, \"D\" -> %" PRID ", \"L\" -> %" PRIL ", \"NV\" -> %" PRIv ", \"NE\" -> %" PRIv ", \"T\" -> %.15f, \"H\" -> {", timestamp, D, L, L * L, D * L * L, T);

  for (q_t i = 0; i < 2; i++) {
    fprintf(outfile_info, "%.15f", H_vec[i]);
    if (i < 2 - 1) {
      fprintf(outfile_info, ", ");
    }
  }

  fprintf(outfile_info, "}, \"GENERATOR\" -> \"%s\", \"EPS\" -> %g |>\n", pert_type, epsilon);

  fclose(outfile_info);

  char *filename_M = (char *)malloc(255 * sizeof(char));
  char *filename_E = (char *)malloc(255 * sizeof(char));
  char *filename_S = (char *)malloc(255 * sizeof(char));
  char *filename_X = (char *)malloc(255 * sizeof(char));

  sprintf(filename_M, "wolff_%lu_M.dat", timestamp);
  sprintf(filename_E, "wolff_%lu_E.dat", timestamp);
  sprintf(filename_S, "wolff_%lu_S.dat", timestamp);
  sprintf(filename_X, "wolff_%lu_X.dat", timestamp);

  FILE *outfile_M = fopen(filename_M, "wb");
  FILE *outfile_E = fopen(filename_E, "wb");
  FILE *outfile_S = fopen(filename_S, "wb");
  FILE *outfile_X = fopen(filename_X, "wb");

  free(filename_M);
  free(filename_E);
  free(filename_S);
  free(filename_X);

  std::function <void(const planar_t *)> *measurements = (std::function <void(const planar_t *)> *)calloc(4, sizeof(std::function <void(const planar_t *)>));

  measurements[0] = (std::function <void(const planar_t *)>)[&](const planar_t *s) {
    float smaller_E = (float)s->E;
    fwrite(&smaller_E, sizeof(float), 1, outfile_E);
  };
  measurements[1] = [&](const planar_t *s) {
    float smaller_X = (float)correlation_length(s);
    fwrite(&smaller_X, sizeof(float), 1, outfile_X);
  };
  measurements[2] = [&](const planar_t *s) {
    write_magnetization(s->M, outfile_M);
  };
  measurements[3] = [&](const planar_t *s) {
    fwrite(&(s->last_cluster_size), sizeof(uint32_t), 1, outfile_S);
  };

  std::function <double(vector_t <2, double>)> H;

  if (modulated_field) {
    H = std::bind(H_modulated, std::placeholders::_1, order, H_vec[0]);
  } else {
    H = std::bind(H_vector <2, double>, std::placeholders::_1, H_vec);
  }

  wolff <orthogonal_t <2, double>, vector_t <2, double>> (N, D, L, T, dot <2, double>, H, gen_R, 4, measurements, silent);

  free(measurements);

  fclose(outfile_M);
  fclose(outfile_E);
  fclose(outfile_S);
  fclose(outfile_X);

  free(H_vec);

  fftw_cleanup();

  return 0;
}

