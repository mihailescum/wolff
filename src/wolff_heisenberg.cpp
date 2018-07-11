
#include <getopt.h>

#include <wolff.h>

int main(int argc, char *argv[]) {

  count_t N = (count_t)1e7;

  D_t D = 2;
  L_t L = 128;
  double T = 2.26918531421;
  double *H = (double *)calloc(MAX_Q, sizeof(double));

  bool silent = false;

  int opt;
  q_t J_ind = 0;
  q_t H_ind = 0;

  while ((opt = getopt(argc, argv, "N:q:D:L:T:J:H:s")) != -1) {
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

  FILE *outfile_info = fopen("wolff_metadata.txt", "a");

  fprintf(outfile_info, "<| \"ID\" -> %lu, \"MODEL\" -> \"HEISENBERG\", q -> \"3\", \"D\" -> %" PRID ", \"L\" -> %" PRIL ", \"NV\" -> %" PRIv ", \"NE\" -> %" PRIv ", \"T\" -> %.15f, \"H\" -> {", timestamp, D, L, L * L, D * L * L, T);

  for (q_t i = 0; i < 3; i++) {
    fprintf(outfile_info, "%.15f", H[i]);
    if (i < 3 - 1) {
      fprintf(outfile_info, ", ");
    }
  }

  fprintf(outfile_info, "} |>\n");

  fclose(outfile_info);


  wolff <orthogonal_t <3, double>, vector_t <3, double>> (N, D, L, T, dot <3, double>, std::bind(H_vector <3, double>, std::placeholders::_1, H), timestamp, silent);

  free(H);

  return 0;
}

