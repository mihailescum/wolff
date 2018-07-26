
#include <types.h>
#include <cmath>
#include <cstring>
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

double squared_mean(int N, double *data) {
  double total = 0;
  for (int i = 0; i < N; i++) {
    total += data[i];
  }

  return total / N;
}

double central_moment(int N, double *data, double mean, int m) {
  double total = 0;
  for (int i = 0; i < N; i++) {
    total += pow(data[i] - mean, m);
  }

  return total / N;
}

void compute_OO(int N, fftw_plan forward_plan, double *forward_data, fftw_plan reverse_plan, double *reverse_data) {

  fftw_execute(forward_plan);

  reverse_data[0] = forward_data[0] * forward_data[0];
  reverse_data[N / 2] = forward_data[N/2] * forward_data[N/2];

  for (count_t i = 1; i < N / 2; i++) {
    reverse_data[i] = pow(forward_data[i], 2) + pow(forward_data[N - i], 2);
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

  fftw_set_timelimit(1);

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

   bool is_finite = 0 == strcmp(model, "ISING") || 0 == strcmp(model, "POTTS") || 0 == strcmp(model, "CLOCK");

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
      } else {
        if (N % 2 == 1 && drop % 2 == 0) {
          drop++; // make sure M is even
        }
      }

      if (N <= drop) {
        printf("\033[F%lu: Number of steps %lu is less than %" PRIcount ", nothing done.\n", id, N, drop);
      } else {
        int M = N - drop;

        double M_f = (double)M;

        if (length > M) {
          length = M;
        }

        double *forward_data = (double *)fftw_malloc(M * sizeof(double));
        fftw_plan forward_plan = fftw_plan_r2r_1d(M, forward_data, forward_data, FFTW_R2HC, 0);
        double *reverse_data = (double *)fftw_malloc(M * sizeof(double));
        fftw_plan reverse_plan = fftw_plan_r2r_1d(M, reverse_data, reverse_data, FFTW_HC2R, 0);


        uint32_t *data_S = (uint32_t *)malloc(N * sizeof(uint32_t));
        fread(data_S, N, sizeof(uint32_t), file_S);
        for (count_t i = 0; i < M; i++) {
          forward_data[i] = (double)data_S[drop + i];
        }
        free(data_S);
        double mean_S = mean(M, forward_data);
        double squaredMean_S = squared_mean(M, forward_data);
        double moment2_S = central_moment(M, forward_data, mean_S, 2);
        double moment4_S = central_moment(M, forward_data, mean_S, 4);

        compute_OO(M, forward_plan, forward_data, reverse_plan, reverse_data);

        sprintf(filename_S, "wolff_%lu_S_OO.dat", id);

        FILE *file_S = fopen(filename_S, "wb");
        fwrite(&M_f, sizeof(double), 1, file_S);
        fwrite(&mean_S, sizeof(double), 1, file_S);
        fwrite(&squaredMean_S, sizeof(double), 1, file_S);
        fwrite(&moment2_S, sizeof(double), 1, file_S);
        fwrite(&moment4_S, sizeof(double), 1, file_S);
        fwrite(reverse_data, sizeof(double), length, file_S);
        fclose(file_S);

        uint32_t *data_B = (uint32_t *)malloc((nb - 1) * N * sizeof(uint32_t));
        uint32_t *data_M = (uint32_t *)malloc((q - 1) * N * sizeof(uint32_t));
        fread(data_B, N * (nb - 1), sizeof(uint32_t), file_B);
        fread(data_M, N * (q - 1), sizeof(uint32_t), file_M);

        for (count_t i = 0; i < M; i++) {
          forward_data[i] = finite_energy(nb, J, q, H, nv, ne, data_B + (nb - 1) * (drop + i), data_M + (q - 1) * (drop + i));
        }

        double mean_E = mean(M, forward_data);
        double squaredMean_E = squared_mean(M, forward_data);
        double moment2_E = central_moment(M, forward_data, mean_E, 2);
        double moment4_E = central_moment(M, forward_data, mean_E, 4);

        free(data_B);
        free(data_M);

        compute_OO(M, forward_plan, forward_data, reverse_plan, reverse_data);

        sprintf(filename_B, "wolff_%lu_E_OO.dat", id);

        FILE *file_E = fopen(filename_B, "wb");
        fwrite(&M_f, sizeof(double), 1, file_E);
        fwrite(&mean_E, sizeof(double), 1, file_E);
        fwrite(&squaredMean_E, sizeof(double), 1, file_E);
        fwrite(&moment2_E, sizeof(double), 1, file_E);
        fwrite(&moment4_E, sizeof(double), 1, file_E);
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
      char *junk = (char *)malloc(1024 * sizeof(char));
      fscanf(metadata, "%[^\n]\n", junk); // throw away the rest of the line, we don't need it
      free(junk);

      char *filename_E = (char *)malloc(128 * sizeof(char));
      char *filename_F = (char *)malloc(128 * sizeof(char));
      char *filename_M = (char *)malloc(128 * sizeof(char));
      char *filename_S = (char *)malloc(128 * sizeof(char));

      sprintf(filename_E, "wolff_%lu_E.dat", id);
      sprintf(filename_F, "wolff_%lu_F.dat", id);
      sprintf(filename_M, "wolff_%lu_M.dat", id);
      sprintf(filename_S, "wolff_%lu_S.dat", id);

      FILE *file_E = fopen(filename_E, "rb");
      FILE *file_F = fopen(filename_F, "rb");
      FILE *file_M = fopen(filename_M, "rb");
      FILE *file_S = fopen(filename_S, "rb");

      fseek(file_S, 0, SEEK_END);
      unsigned long N = ftell(file_S) / sizeof(uint32_t);
      fseek(file_S, 0, SEEK_SET);

      if (speedy_drop) {
        drop = N - pow(2, floor(log(N) / log(2)));
      } else {
        if (N % 2 == 1 && drop % 2 == 0) {
          drop++; // make sure M is even
        }
      }

      if (N <= drop) {
        printf("\033[F%lu: Number of steps %lu is less than %" PRIcount ", nothing done.\n", id, N, drop);
      } else {
        int M = N - drop;
        double M_f = (double)M;

        if (length > M) {
          length = M;
        }

        double *forward_data = (double *)fftw_malloc(M * sizeof(double));
        fftw_plan forward_plan = fftw_plan_r2r_1d(M, forward_data, forward_data, FFTW_R2HC, 0);

        double *reverse_data = (double *)fftw_malloc(M * sizeof(double));
        fftw_plan reverse_plan = fftw_plan_r2r_1d(M, reverse_data, reverse_data, FFTW_HC2R, 0);

        if (file_S != NULL) {
          uint32_t *data_S = (uint32_t *)malloc(N * sizeof(uint32_t));

          fread(data_S, sizeof(uint32_t), N, file_S);
          fclose(file_S);

          for (int i = 0; i < M; i++) {
            forward_data[i] = (double)data_S[drop + i];
          }
          free(data_S);

          double mean_S = mean(M, forward_data);
          double squaredMean_S = squared_mean(M, forward_data);
          double moment2_S = central_moment(M, forward_data, mean_S, 2);
          double moment4_S = central_moment(M, forward_data, mean_S, 4);

          compute_OO(M, forward_plan, forward_data, reverse_plan, reverse_data);

          sprintf(filename_S, "wolff_%lu_S_OO.dat", id);
          FILE *file_S_new = fopen(filename_S, "wb");
          fwrite(&M_f, sizeof(double), 1, file_S_new);
          fwrite(&mean_S, sizeof(double), 1, file_S_new);
          fwrite(&squaredMean_S, sizeof(double), 1, file_S_new);
          fwrite(&moment2_S, sizeof(double), 1, file_S_new);
          fwrite(&moment4_S, sizeof(double), 1, file_S_new);
          fwrite(reverse_data, sizeof(double), length, file_S_new);
          fclose(file_S_new);
        }
        if (file_F != NULL) {
          float *data_F = (float *)malloc(N * sizeof(float));

          fread(data_F, sizeof(float), N, file_F);
          fclose(file_F);

          for (int i = 0; i < M; i++) {
            forward_data[i] = (double)data_F[drop + i];
          }
          free(data_F);

          double mean_F = mean(M, forward_data);
          double squaredMean_F = squared_mean(M, forward_data);
          double moment2_F = central_moment(M, forward_data, mean_F, 2);
          double moment4_F = central_moment(M, forward_data, mean_F, 4);

          compute_OO(M, forward_plan, forward_data, reverse_plan, reverse_data);

          sprintf(filename_F, "wolff_%lu_F_OO.dat", id);
          FILE *file_F_new = fopen(filename_F, "wb");
          fwrite(&M_f, sizeof(double), 1, file_F_new);
          fwrite(&mean_F, sizeof(double), 1, file_F_new);
          fwrite(&squaredMean_F, sizeof(double), 1, file_F_new);
          fwrite(&moment2_F, sizeof(double), 1, file_F_new);
          fwrite(&moment4_F, sizeof(double), 1, file_F_new);
          fwrite(reverse_data, sizeof(double), length, file_F_new);
          fclose(file_F_new);
        }
        if (file_E != NULL) {
          float *data_E = (float *)malloc(N * sizeof(float));

          fread(data_E, sizeof(float), N, file_E);
          fclose(file_E);

          for (int i = 0; i < M; i++) {
            forward_data[i] = (double)data_E[drop + i];
          }
          free(data_E);

          double mean_E = mean(M, forward_data);
          double squaredMean_E = squared_mean(M, forward_data);
          double moment2_E = central_moment(M, forward_data, mean_E, 2);
          double moment4_E = central_moment(M, forward_data, mean_E, 4);

          compute_OO(M, forward_plan, forward_data, reverse_plan, reverse_data);

          sprintf(filename_E, "wolff_%lu_E_OO.dat", id);
          FILE *file_E_new = fopen(filename_E, "wb");
          fwrite(&M_f, sizeof(double), 1, file_E_new);
          fwrite(&mean_E, sizeof(double), 1, file_E_new);
          fwrite(&squaredMean_E, sizeof(double), 1, file_E_new);
          fwrite(&moment2_E, sizeof(double), 1, file_E_new);
          fwrite(&moment4_E, sizeof(double), 1, file_E_new);
          fwrite(reverse_data, sizeof(double), length, file_E_new);
          fclose(file_E_new);
        }
        if (file_M != NULL) {
          if (0 == strcmp(model, "PLANAR")) {
            float *data_M = (float *)malloc(2 * N * sizeof(float));
            fread(data_M, sizeof(float), 2 * N, file_M);
            fclose(file_M);
            for (int i = 0; i < M; i++) {
              forward_data[i] = (double)sqrt(pow(data_M[2 * drop + 2 * i], 2) + pow(data_M[2 * drop + 2 * i + 1], 2));
            }
            free(data_M);
          } else if (0 == strcmp(model, "HEISENBERG")) {
            float *data_M = (float *)malloc(3 * N * sizeof(float));
            fread(data_M, sizeof(float), 3 * N, file_M);
            fclose(file_M);
            for (int i = 0; i < M; i++) {
              forward_data[i] = sqrt(pow(data_M[3 * drop + 3 * i], 2) + pow(data_M[3 * drop + 3 * i + 1], 2) + pow(data_M[3 * drop + 3 * i + 2], 2));
            }
            free(data_M);
          } else {
            printf("UNKNOWN MODEL\n");
            exit(EXIT_FAILURE);
          }

          double mean_M = mean(M, forward_data);
          double squaredMean_M = squared_mean(M, forward_data);
          double moment2_M = central_moment(M, forward_data, mean_M, 2);
          double moment4_M = central_moment(M, forward_data, mean_M, 4);

          compute_OO(M, forward_plan, forward_data, reverse_plan, reverse_data);

          sprintf(filename_M, "wolff_%lu_M_OO.dat", id);
          FILE *file_M_new = fopen(filename_M, "wb");
          fwrite(&M_f, sizeof(double), 1, file_M_new);
          fwrite(&mean_M, sizeof(double), 1, file_M_new);
          fwrite(&squaredMean_M, sizeof(double), 1, file_M_new);
          fwrite(&moment2_M, sizeof(double), 1, file_M_new);
          fwrite(&moment4_M, sizeof(double), 1, file_M_new);
          fwrite(reverse_data, sizeof(double), length, file_M_new);
          fclose(file_M_new);
        }

        printf("\033[F%lu: Correlation functions for %d steps written.\n", id, M);
        fftw_destroy_plan(forward_plan);
        fftw_destroy_plan(reverse_plan);
        fftw_free(forward_data);
        fftw_free(reverse_data);

      }
      free(filename_E);
      free(filename_S);
      free(filename_F);
      free(filename_M);
    }
  }

  free(model);
  fclose(metadata);
  fftw_cleanup();

  return 0;
}

