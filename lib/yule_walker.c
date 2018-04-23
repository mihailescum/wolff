
#include "yule_walker.h"

double yule_walker(const autocorr_t *a) {
  gsl_vector *rhos = gsl_vector_alloc(a->W);
  gsl_matrix *mat = gsl_matrix_alloc(a->W, a->W);
  gsl_vector *phis = gsl_vector_alloc(a->W);

  for (count_t i = 0; i < a->W; i++) {
    gsl_vector_set(rhos, i, rho(a, i));
  }

  for (count_t i = 0; i < a->W; i++) {
    gsl_matrix_set(mat, i, i, 1.0);

    for (count_t j = 0; j < i; j++) {
      gsl_matrix_set(mat, i, j, gsl_vector_get(rhos, i - 1 - j));
    }

    for (count_t j = 0; j < a->W - 1 - i; j++) {
      gsl_matrix_set(mat, i, i + 1 + j, gsl_vector_get(rhos, j));
    }
  }

  int out = gsl_linalg_cholesky_solve(mat, rhos, phis);

  double rhophi = 0;
  double onephi = 0;

  for (count_t i = 0; i < a->W; i++) {
    rhophi += gsl_vector_get(rhos, i) * gsl_vector_get(phis, i);
    onephi += gsl_vector_get(phis, i);
  }

  gsl_vector_free(rhos);
  gsl_matrix_free(mat);
  gsl_vector_free(phis);

  return (1 - rhophi) / pow(1 - onephi, 2);
}

