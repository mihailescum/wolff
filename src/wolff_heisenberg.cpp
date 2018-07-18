
#include <getopt.h>

#include <correlation.h>
#include <wolff.h>

typedef state_t <orthogonal_t <3, double>, vector_t <3, double>> sim_t;

int main(int argc, char *argv[]) {

  count_t N = (count_t)1e7;

  D_t D = 2;
  L_t L = 128;
  double T = 2.26918531421;
  double *H = (double *)calloc(MAX_Q, sizeof(double));

  bool silent = false;
  bool use_pert = false;

  int opt;
  q_t J_ind = 0;
  q_t H_ind = 0;
  double epsilon = 1;

  while ((opt = getopt(argc, argv, "N:q:D:L:T:J:H:spe:")) != -1) {
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
      H[H_ind] = atof(optarg);
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

  std::function <orthogonal_t <3, double>(gsl_rng *, const sim_t *)> gen_R;

  const char *pert_type;

  if (use_pert) {
    gen_R = std::bind(generate_rotation_perturbation <3>, std::placeholders::_1, std::placeholders::_2, epsilon);
    pert_type = "PERTURB";
  } else {
    gen_R = generate_rotation_uniform <3>;
    pert_type = "UNIFORM";
  }

  FILE *outfile_info = fopen("wolff_metadata.txt", "a");

  fprintf(outfile_info, "<| \"ID\" -> %lu, \"MODEL\" -> \"HEISENBERG\", \"q\" -> 3, \"D\" -> %" PRID ", \"L\" -> %" PRIL ", \"NV\" -> %" PRIv ", \"NE\" -> %" PRIv ", \"T\" -> %.15f, \"H\" -> {", timestamp, D, L, L * L, D * L * L, T);

  for (q_t i = 0; i < 3; i++) {
    fprintf(outfile_info, "%.15f", H[i]);
    if (i < 3 - 1) {
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

  std::function <void(const sim_t *)> *measurements = (std::function <void(const sim_t *)> *)malloc(4 * sizeof(std::function <void(const sim_t *)>));

  measurements[0] = [&](const sim_t *s) {
    float smaller_E = (float)s->E;
    fwrite(&smaller_E, sizeof(float), 1, outfile_E);
  };
  measurements[1] = [&](const sim_t *s) {
    float smaller_X = (float)correlation_length(s);
    fwrite(&smaller_X, sizeof(float), 1, outfile_X);
  };
  measurements[2] = [&](const sim_t *s) {
    write_magnetization(s->M, outfile_M);
  };
  measurements[3] = [&](const sim_t *s) {
    fwrite(&(s->last_cluster_size), sizeof(uint32_t), 1, outfile_S);
  };

  wolff <orthogonal_t <3, double>, vector_t <3, double>> (N, D, L, T, dot <3, double>, std::bind(H_vector <3, double>, std::placeholders::_1, H), gen_R, 4, measurements, silent);

  free(measurements);

  fclose(outfile_M);
  fclose(outfile_E);
  fclose(outfile_S);
  fclose(outfile_X);

  free(H);

  fftw_cleanup();

  return 0;
}

