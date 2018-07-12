
#include <types.h>
#include <measurement.h>
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <fftw3.h>

template <class T>
double mean(int N, T *data) {
  double total = 0;
  for (int i = 0; i < N; i++) {
    total += (double)data[i];
  }

  return total / N;
}

void compute_OO(int N, fftw_plan forward_plan, double *forward_data, fftw_plan reverse_plan, double *reverse_data) {

  fftw_execute(forward_plan);

  reverse_data[0] = forward_data[0] * forward_data[0];
  reverse_data[N / 2] = forward_data[N/2] * forward_data[N/2];

  for (count_t i = 1; i < N / 2 - 1; i++) {
    reverse_data[i] = forward_data[i] * forward_data[i] + forward_data[N - i] * forward_data[N - i];
    reverse_data[N - i] = 0;
  }

  fftw_execute(reverse_plan);
}

double finite_energy(q_t nb, double *J, q_t q, double *H, v_t nv, v_t ne, uint32_t *bo, uint32_t *so) {
  double energy = 0;

  v_t tot = 0;
  for (q_t i = 0; i < nb - 1; i++) {
    energy -= J[i] * bo[i];
    tot += bo[i];
  }

  energy -= J[nb - 1] * (ne - tot);

  tot = 0;
  for (q_t i = 0; i < q - 1; i++) {
    energy -= H[i] * so[i];
    tot += so[i];
  }

  energy -= H[q - 1] * (nv - tot);

  return energy;
}

int main (int argc, char *argv[]) {
  count_t drop = (count_t)1e4;
  count_t length = (count_t)1e4;
  bool speedy_drop = false;
  bool from_stdin = false;

  int opt;

  while ((opt = getopt(argc, argv, "d:l:sp")) != -1) {
    switch (opt) {
      case 'd':
        drop = (count_t)atof(optarg);
        break;
      case 'l':
        length = (count_t)atof(optarg);
        break;
      case 's':
        speedy_drop = true;
        break;
      case 'p':
        from_stdin = true;
        break;
      default:
        exit(EXIT_FAILURE);
    }
  }
  FILE *metadata;

  fftw_set_timelimit(10);

  if (from_stdin) {
    metadata = stdin;
  } else {
    metadata = fopen("wolff_metadata.txt", "r");
  }

  if (metadata == NULL) {
    printf("Metadata file not found. Make sure you are in the correct directory!\n");
    exit(EXIT_FAILURE);
  }

  unsigned long id;
  char *model = (char *)malloc(32 * sizeof(char));

  if (model == NULL) {
    printf("Malloc failed.\n");
    exit(EXIT_FAILURE);
  }

  q_t q;
  D_t D;
  L_t L;
  v_t nv, ne;

 while (EOF != fscanf(metadata, "<| \"ID\" -> %lu, \"MODEL\" -> \"%[^\"]\", \"q\" -> %" SCNq ", \"D\" -> %" SCND ", \"L\" -> %" SCNL ", \"NV\" -> %" SCNv ", \"NE\" -> %" SCNv ", ", &id, model, &q, &D, &L, &nv, &ne)) {

   printf("%lu: Processing...\n", id);

   bool is_finite = 0 == strcmp(model, "ISING") || 0 == strcmp(model, "POTTS") || 0 == strcmp(model, "CLOCK") || 0 == strcmp(model, "DGM");

    if (is_finite) {
      q_t nb;
      double T;
      fscanf(metadata, "\"NB\" -> %" SCNq ", \"T\" -> %lf, \"J\" -> {", &nb, &T);
      double *J = (double *)malloc(nb * sizeof(double));
      double *H = (double *)malloc(q * sizeof(double));

      if (J == NULL || H == NULL) {
        printf("%lu: Malloc failed.\n", id);
        break;
      }

      for (q_t i = 0; i < nb - 1; i++) {
        fscanf(metadata, "%lf, ", &(J[i]));
      }
      fscanf(metadata, "%lf}, \"H\" -> {", &(J[nb - 1]));
      for (q_t i = 0; i < q - 1; i++) {
        fscanf(metadata, "%lf, ", &(H[i]));
      }
      fscanf(metadata, "%lf} |>\n", &(H[q - 1]));

      char *filename_M = (char *)malloc(128 * sizeof(char));
      char *filename_B = (char *)malloc(128 * sizeof(char));
      char *filename_S = (char *)malloc(128 * sizeof(char));

      if (filename_M == NULL || filename_B == NULL || filename_S == NULL) {
        printf("%lu: Malloc failed.\n", id);
        break;
      }

      sprintf(filename_M, "wolff_%lu_M.dat", id);
      sprintf(filename_B, "wolff_%lu_B.dat", id);
      sprintf(filename_S, "wolff_%lu_S.dat", id);

      FILE *file_M = fopen(filename_M, "rb");
      FILE *file_B = fopen(filename_B, "rb");
      FILE *file_S = fopen(filename_S, "rb");

      if (file_M == NULL || file_B == NULL || file_S == NULL) {
        printf("%lu: Opening data file failed.\n", id);
        break;
      }

      fseek(file_S, 0, SEEK_END);
      unsigned long N = ftell(file_S) / sizeof(uint32_t);
      fseek(file_S, 0, SEEK_SET);

      if (speedy_drop) {
        drop = N - pow(2, floor(log(N) / log(2)));
      }

      if (N < drop) {
        printf("%lu: Number of steps %lu is less than %" PRIcount ", nothing done.\n", id, N, drop);
        break;
      } else {

        uint32_t *data_S = (uint32_t *)malloc(N * sizeof(uint32_t));
        uint32_t *data_B = (uint32_t *)malloc((nb - 1) * N * sizeof(uint32_t));
        uint32_t *data_M = (uint32_t *)malloc((q - 1) * N * sizeof(uint32_t));

        fread(data_S, N, sizeof(uint32_t), file_S);
        fread(data_B, N * (nb - 1), sizeof(uint32_t), file_B);
        fread(data_M, N * (q - 1), sizeof(uint32_t), file_M);

        int M = N - drop;

        double *forward_data = (double *)fftw_malloc(M * sizeof(double));
        fftw_plan forward_plan = fftw_plan_r2r_1d(M, forward_data, forward_data, FFTW_R2HC, 0);
        double *reverse_data = (double *)fftw_malloc(M * sizeof(double));
        fftw_plan reverse_plan = fftw_plan_r2r_1d(M, reverse_data, reverse_data, FFTW_HC2R, 0);

        double mean_S = mean(M, data_S);
        double M_f = (double)M;

        for (count_t i = 0; i < M; i++) {
          forward_data[i] = finite_energy(nb, J, q, H, nv, ne, data_B + (nb - 1) * (drop + i), data_M + (q - 1) * (drop + i));
        }

        double mean_E = mean(M, forward_data);

        free(data_B);
        free(data_M);
        free(data_S);

        compute_OO(M, forward_plan, forward_data, reverse_plan, reverse_data);

        sprintf(filename_B, "wolff_%lu_E_OO.dat", id);

        FILE *file_E = fopen(filename_B, "wb");
        fwrite(&M_f, sizeof(double), 1, file_E);
        fwrite(&mean_E, sizeof(double), 1, file_E);
        fwrite(&mean_S, sizeof(double), 1, file_E);
        fwrite(reverse_data, sizeof(double), length, file_E);
        fclose(file_E);

        printf("\033[F%lu: Correlation functions for %d steps written.\n", id, M);

        fftw_destroy_plan(forward_plan);
        fftw_destroy_plan(reverse_plan);
        fftw_free(forward_data);
        fftw_free(reverse_data);

      }

      fclose(file_M);
      fclose(file_B);
      fclose(file_S);

      free(J);
      free(H);

      free(filename_S);
      free(filename_B);
      free(filename_M);

    } else {
      double T;
      fscanf(metadata, "\"T\" -> %lf, \"H\" -> {", &T);
      double *H = (double *)malloc(q * sizeof(double));

      for (q_t i = 0; i < q - 1; i++) {
        fscanf(metadata, "%lf, ", &(H[i]));
      }
      fscanf(metadata, "%lf} |>\n", &(H[q - 1]));
      free(H);

      char *filename_E = (char *)malloc(128 * sizeof(char));
      char *filename_S = (char *)malloc(128 * sizeof(char));

      sprintf(filename_E, "wolff_%lu_E.dat", id);
      sprintf(filename_S, "wolff_%lu_S.dat", id);

      FILE *file_E = fopen(filename_E, "rb");
      FILE *file_S = fopen(filename_S, "rb");

      fseek(file_S, 0, SEEK_END);
      unsigned long N = ftell(file_S) / sizeof(uint32_t);
      fseek(file_S, 0, SEEK_SET);

      if (speedy_drop) {
        drop = N - pow(2, floor(log(N) / log(2)));
      }

      if (N < drop) {
        printf("%lu: Number of steps %lu is less than %" PRIcount ", nothing done.\n", id, N, drop);
      } else {
        int M = N - drop;
        double M_f = (double)M;

        uint32_t *data_S = (uint32_t *)malloc(N * sizeof(uint32_t));

        fread(data_S, sizeof(uint32_t), N, file_S);

        double *forward_data = (double *)fftw_malloc(M * sizeof(double));
        double *reverse_data = (double *)fftw_malloc(M * sizeof(double));

        fftw_plan forward_plan = fftw_plan_r2r_1d(M, forward_data, forward_data, FFTW_R2HC, 0);
        fftw_plan reverse_plan = fftw_plan_r2r_1d(M, reverse_data, reverse_data, FFTW_HC2R, 0);

        double mean_S = mean(M, data_S);

        free(data_S);

        float *data_E = (float *)malloc(N * sizeof(float));
        fread(data_E, sizeof(float), N, file_E);
        for (int i = 0; i < M; i++) {
          forward_data[i] = (double)data_E[drop + i];
        }

        double mean_E = mean(M, forward_data);

        compute_OO(M, forward_plan, forward_data, reverse_plan, reverse_data);

        sprintf(filename_E, "wolff_%lu_E_OO.dat", id);
        FILE *file_E_new = fopen(filename_E, "wb");
        fwrite(&M_f, sizeof(double), 1, file_E_new);
        fwrite(&mean_E, sizeof(double), 1, file_E_new);
        fwrite(&mean_S, sizeof(double), 1, file_E_new);
        fwrite(reverse_data, sizeof(double), length, file_E_new);
        fclose(file_E_new);

        free(data_E);
        printf("\033[F%lu: Correlation functions for %d steps written.\n", id, M);
        fftw_destroy_plan(forward_plan);
        fftw_destroy_plan(reverse_plan);
        fftw_free(forward_data);
        fftw_free(reverse_data);

      }
      free(filename_E);
      free(filename_S);

      fclose(file_E);
      fclose(file_S);

    }
  }

  free(model);
  fclose(metadata);
  fftw_cleanup();

  return 0;
}

