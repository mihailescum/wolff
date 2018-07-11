
#include <types.h>
#include <measurement.h>
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>

template <class T>
double *compute_OO(count_t length, T *data, count_t N) {
  double *OO = (double *)calloc(2 + length, sizeof(double));
  OO[0] = (double)N;

  for (count_t i = 0; i < N; i++) {
    OO[1] += data[i];
    for (count_t j = 0; j < length; j++) {
      if ((long int)N - (long int)i - (long int)j > 0) {
        OO[2 + j] += data[i] * data[i + j]; 
      }
    }
  }

  OO[1] /= N;

  for (count_t j = 0; j < length; j++) {
    OO[2 + j] /= N - j;
  }

  return OO;
}

float finite_energy(q_t nb, double *J, q_t q, double *H, v_t nv, v_t ne, uint32_t *bo, uint32_t *so) {
  float energy = 0;

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

        uint32_t *data_S = (uint32_t *)malloc(N * sizeof(uint32_t));
        uint32_t *data_B = (uint32_t *)malloc((nb - 1) * N * sizeof(uint32_t));
        uint32_t *data_M = (uint32_t *)malloc((q - 1) * N * sizeof(uint32_t));

        fread(data_S, N, sizeof(uint32_t), file_S);
        fread(data_B, N * (nb - 1), sizeof(uint32_t), file_B);
        fread(data_M, N * (q - 1), sizeof(uint32_t), file_M);

        double *OO_S = compute_OO(length, data_S + drop, N - drop);

        sprintf(filename_S, "wolff_%lu_S_OO.dat", id);

        FILE *file_S_new = fopen(filename_S, "wb");
        fwrite(OO_S, 2 + length, sizeof(double), file_S_new);
        fclose(file_S_new);
        free(data_S);

        float *data_E = (float *)malloc(N * sizeof(float));

        for (count_t i = 0; i < N; i++) {
          data_E[i] = finite_energy(nb, J, q, H, nv, ne, data_B + (nb - 1) * i, data_M + (q - 1) * i);
        }

        free(data_B);
        free(data_M);

        double *OO_E = compute_OO(length, data_E + drop, N - drop);

        free(data_E);
        sprintf(filename_B, "wolff_%lu_E_OO.dat", id);

        FILE *file_E = fopen(filename_B, "wb");
        fwrite(OO_E, 2 + length, sizeof(double), file_E);
        fclose(file_E);

        free(OO_S);
        free(OO_E);

        printf("Correlation functions for %lu steps written.\n", N - drop);
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

        uint32_t *data_S = (uint32_t *)malloc(N * sizeof(uint32_t));
        float *data_E = (float *)malloc(N * sizeof(float));

        fread(data_S, N, sizeof(uint32_t), file_S);
        fread(data_E, N, sizeof(float), file_E);

        double *OO_S = compute_OO(length, data_S + drop, N - drop);
        double *OO_E = compute_OO(length, data_E + drop, N - drop);

        sprintf(filename_S, "wolff_%lu_S_OO.dat", id);
        sprintf(filename_E, "wolff_%lu_E_OO.dat", id);

        FILE *file_S_new = fopen(filename_S, "wb");
        FILE *file_E_new = fopen(filename_E, "wb");
        fwrite(OO_S, 1 + length, sizeof(double), file_S_new);
        fwrite(OO_E, 1 + length, sizeof(double), file_E_new);
        fclose(file_S_new);
        fclose(file_E_new);
        free(data_S);
        free(data_E);
        free(OO_S);
        free(OO_E);

        printf("Correlation functions for %lu steps written.\n", N - drop);
      }
      free(filename_E);
      free(filename_S);

      fclose(file_E);
      fclose(file_S);

    }
  }
  free(model);
  fclose(metadata);

  return 0;
}

