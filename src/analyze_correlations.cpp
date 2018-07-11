
#include <types.h>
#include <measurement.h>
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <fftw3.h>

double *compute_OO(count_t length, double *tmp_data, count_t N) {
  double *OO = (double *)calloc(2 + length, sizeof(double));

  if (OO == NULL) {
    printf("calloc failure\n");
    exit(EXIT_FAILURE);
  }

//  int better_N = pow(2, floor(log(N) / log(2)));
  int better_N = N;
  double *data = tmp_data + (N - better_N);

  OO[0] = (double)better_N;

  for (count_t i = 0; i < better_N; i++) {
    OO[1] += data[i];
  }

  OO[1] /= better_N;

  fftw_plan plan = fftw_plan_r2r_1d(better_N, data, data, FFTW_R2HC, FFTW_ESTIMATE);
  fftw_execute(plan);

  double *tmp = (double *)malloc(better_N * sizeof(double));

  tmp[0] = data[0] * data[0];
  tmp[better_N / 2] = data[better_N/2] * data[better_N/2];

  for (count_t i = 1; i < better_N / 2 - 1; i++) {
    tmp[i] = data[i] * data[i] + data[better_N - i] * data[better_N - i];
    tmp[better_N - i] = 0;
  }

  fftw_destroy_plan(plan);

  fftw_plan plan2 = fftw_plan_r2r_1d(better_N, tmp, tmp, FFTW_HC2R, FFTW_ESTIMATE);
  fftw_execute(plan2);

  for (count_t j = 0; j < length; j++) {
    OO[2 + j] = tmp[j] / pow(better_N, 2);
  }

  free(tmp);
  fftw_destroy_plan(plan2);

  return OO;
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

  int opt;

  while ((opt = getopt(argc, argv, "d:l:")) != -1) {
    switch (opt) {
      case 'd':
        drop = (count_t)atof(optarg);
        break;
      case 'l':
        length = (count_t)atof(optarg);
        break;
      default:
        exit(EXIT_FAILURE);
    }
  }

  FILE *metadata = fopen("wolff_metadata.txt", "r");

  if (metadata == NULL) {
    exit(EXIT_FAILURE);
  }

  unsigned long id;
  char *model = (char *)malloc(32 * sizeof(char));
  q_t q;
  D_t D;
  L_t L;
  v_t nv, ne;

 while (EOF != fscanf(metadata, "<| \"ID\" -> %lu, \"MODEL\" -> \"%[^\"]\", \"q\" -> %" SCNq ", \"D\" -> %" SCND ", \"L\" -> %" SCNL ", \"NV\" -> %" SCNv ", \"NE\" -> %" SCNv ", ", &id, model, &q, &D, &L, &nv, &ne)) {

    if (0 == strcmp(model, "ISING") || 0 == strcmp(model, "POTTS") || 0 == strcmp(model, "CLOCK") || 0 == strcmp(model, "DGM")) {
      q_t nb;
      double T;
      fscanf(metadata, "\"NB\" -> %" SCNq ", \"T\" -> %lf, \"J\" -> {", &nb, &T);
      double *J = (double *)malloc(nb * sizeof(double));
      double *H = (double *)malloc(q * sizeof(double));

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

      sprintf(filename_M, "wolff_%lu_M.dat", id);
      sprintf(filename_B, "wolff_%lu_B.dat", id);
      sprintf(filename_S, "wolff_%lu_S.dat", id);

      FILE *file_M = fopen(filename_M, "rb");
      FILE *file_B = fopen(filename_B, "rb");
      FILE *file_S = fopen(filename_S, "rb");

      fseek(file_S, 0, SEEK_END);
      unsigned long N = ftell(file_S) / sizeof(uint32_t);
      fseek(file_S, 0, SEEK_SET);

      if (N < drop) {
        printf("Number of steps %lu is less than %" PRIcount ", nothing done.\n", N, drop);
      } else {

        double *data_S = (double *)malloc(N * sizeof(double));
        uint32_t *data_B = (uint32_t *)malloc((nb - 1) * N * sizeof(uint32_t));
        uint32_t *data_M = (uint32_t *)malloc((q - 1) * N * sizeof(uint32_t));

        uint32_t tmp;
        for (count_t i = 0; i < N; i++) {
          fread(&tmp, 1, sizeof(uint32_t), file_S);
          data_S[i] = (double)tmp;
        }

        fread(data_B, N * (nb - 1), sizeof(uint32_t), file_B);
        fread(data_M, N * (q - 1), sizeof(uint32_t), file_M);

        double *OO_S = compute_OO(length, data_S + drop, N - drop);

        sprintf(filename_S, "wolff_%lu_S_OO.dat", id);

        FILE *file_S_new = fopen(filename_S, "wb");
        fwrite(OO_S, sizeof(double), 2 + length, file_S_new);
        fclose(file_S_new);
        free(data_S);

        double *data_E = (double *)malloc(N * sizeof(double));

        for (count_t i = 0; i < N; i++) {
          data_E[i] = finite_energy(nb, J, q, H, nv, ne, data_B + (nb - 1) * i, data_M + (q - 1) * i);
        }

        free(data_B);
        free(data_M);

        double *OO_E = compute_OO(length, data_E + drop, N - drop);

        free(data_E);
        sprintf(filename_B, "wolff_%lu_E_OO.dat", id);

        FILE *file_E = fopen(filename_B, "wb");
        if (file_E == NULL) {
          printf("failed to open file\n");
          exit(EXIT_FAILURE);
        }
        fwrite(OO_E, sizeof(double), 2 + length, file_E);
        fclose(file_E);

        printf("Correlation functions for %g steps written.\n", OO_S[0]);

        free(OO_S);
        free(OO_E);

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

      if (N < drop) {
        printf("Number of steps %lu is less than %" PRIcount ", nothing done.\n", N, drop);
      } else {

        double *data_S = (double *)malloc(N * sizeof(double));
        double *data_E = (double *)malloc(N * sizeof(double));

        uint32_t tmp;
        float tmp2;
        for (count_t i = 0; i < N; i++) {
          fread(&tmp, 1, sizeof(uint32_t), file_S);
          data_S[i] = (double)tmp;
          fread(&tmp2, 1, sizeof(float), file_E);
          data_E[i] = (double)tmp2;
        }

        double *OO_S = compute_OO(length, data_S + drop, N - drop);
        double *OO_E = compute_OO(length, data_E + drop, N - drop);

        sprintf(filename_S, "wolff_%lu_S_OO.dat", id);
        sprintf(filename_E, "wolff_%lu_E_OO.dat", id);

        FILE *file_S_new = fopen(filename_S, "wb");
        fwrite(OO_S, sizeof(double), 2 + length, file_S_new);
        fclose(file_S_new);

        FILE *file_E_new = fopen(filename_E, "wb");
        fwrite(OO_E, sizeof(double), 2 + length, file_E_new);
        fclose(file_E_new);

        free(data_S);
        free(data_E);
        printf("Correlation functions for %g steps written.\n", OO_S[0]);
        free(OO_S);
        free(OO_E);

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

