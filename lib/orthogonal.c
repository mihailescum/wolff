
#include "orthogonal.h"

void vector_replace(q_t n, double *v1, const double *v2) {
  for (q_t i = 0; i < n; i++) {
    v1[i] = v2[i];
  }
}

void vector_add(q_t n, double *v1, const double *v2) {
  for (q_t i = 0; i < n; i++) {
    v1[i] += v2[i];
  }
}

void vector_subtract(q_t n, double *v1, const double *v2) {
  for (q_t i = 0; i < n; i++) {
    v1[i] -= v2[i];
  }
}

double *vector_rotate(q_t n, double *rot, double *vec) {
  double *rot_vec = (double *)malloc(n * sizeof(double));

  double prod = 0.0;
  for (q_t i = 0; i < n; i++) {
    prod += rot[i] * vec[i];
  }

  for (q_t i = 0; i < n; i++) {
    rot_vec[i] = vec[i] - 2 * prod * rot[i];
  }

  return rot_vec;
}

double *vector_rotate_inverse(q_t n, const double *rot, const double *vec) {
  double *rot_vec = (double *)calloc(n, sizeof(double));

  for (q_t i = 0; i < n; i++) {
    for (q_t j = 0; j < n; j++) {
      rot_vec[i] += rot[n * j + i] * vec[j];
    }
  }

  return rot_vec;
}

double vector_dot(q_t n, double *v1, double *v2) {
  double dot = 0;

  for (q_t i = 0; i < n; i++) {
    dot += v1[i] * v2[i];
  }

  return dot;
}

double *orthogonal_rotate(q_t n, double *r, double *m) {
  double *mul = (double *)calloc(n * n, sizeof(double));

  for (q_t i = 0; i < n; i++) {
    double akOki = 0;

    for (q_t k = 0; k < n; k++) {
      akOki += r[k] * m[n * k + i];
    }

    for (q_t j = 0; j < n; j++) {
      mul[n * j + i] = m[n * j + i] - 2 * akOki * r[j];
    }
  }

  return mul;
}

double *gen_rot(gsl_rng *r, q_t n) {
  double *v = (double *)malloc(n * sizeof(double));

  double v2 = 0;

  for (q_t i = 0; i < n; i++) {
    v[i] = gsl_ran_ugaussian(r);
    v2 += v[i] * v[i];
  }

  double magv = sqrt(v2);

  for (q_t i = 0; i < n; i++) {
    v[i] /= magv;
  }

  return v;
}

