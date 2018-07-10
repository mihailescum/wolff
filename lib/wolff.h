
#include <time.h>
#include <getopt.h>

#include <cluster.h>

template <q_t q, class T>
double H_vector(vector_t <q, T> v1, T *H) {
  vector_t <q, T> H_vec;
  H_vec.x = H;
  return (double)(dot <q, T> (v1, H_vec));
}

template <class R_t, class X_t>
void wolff(count_t N, D_t D, L_t L, double T, std::function <double(X_t, X_t)> J, std::function <double(X_t)> H, unsigned long timestamp, bool silent) {

  state_t <R_t, X_t> s(D, L, T, J, H);

  // initialize random number generator
  gsl_rng *r = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(r, rand_seed());

  char *filename_M = (char *)malloc(255 * sizeof(char));
  char *filename_E = (char *)malloc(255 * sizeof(char));
  char *filename_S = (char *)malloc(255 * sizeof(char));

  sprintf(filename_M, "wolff_%lu_M.dat", timestamp);
  sprintf(filename_E, "wolff_%lu_E.dat", timestamp);
  sprintf(filename_S, "wolff_%lu_S.dat", timestamp);

  FILE *outfile_M = fopen(filename_M, "wb");
  FILE *outfile_E = fopen(filename_E, "wb");
  FILE *outfile_S = fopen(filename_S, "wb");

  free(filename_M);
  free(filename_E);
  free(filename_S);

  v_t cluster_size = 0;

  if (!silent) printf("\n");
  for (count_t steps = 0; steps < N; steps++) {
    if (!silent) printf("\033[F\033[JWOLFF: sweep %" PRIu64 " / %" PRIu64 ": E = %.2f, S = %" PRIv "\n", steps, N, s.E, cluster_size);

    v_t v0 = gsl_rng_uniform_int(r, s.nv);
     
    R_t step;
    generate_rotation(r, &step);

    cluster_size = flip_cluster <R_t, X_t> (&s, v0, step, r);

    free_spin(step);

    {
      float smaller_E = (float)s.E;
      fwrite(&smaller_E, sizeof(float), 1, outfile_E);
    }
    write_magnetization(s.M, outfile_M);
    fwrite(&cluster_size, sizeof(uint32_t), 1, outfile_S);

  }
  if (!silent) {
    printf("\033[F\033[J");
  }
  printf("WOLFF: sweep %" PRIu64 " / %" PRIu64 ": E = %.2f, S = %" PRIv "\n", N, N, s.E, cluster_size);

  fclose(outfile_M);
  fclose(outfile_E);
  fclose(outfile_S);

  gsl_rng_free(r);
}

