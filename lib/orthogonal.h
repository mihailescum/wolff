
#include <stdlib.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <math.h>

#include "types.h"

void vector_replace(q_t n, double *v1, const double *v2);

void vector_add(q_t n, double *v1, const double *v2);

void vector_subtract(q_t n, double *v1, const double *v2);

double *vector_rotate(q_t n, const double *rot, const double *vec);

double *vector_rotate_inverse(q_t n, const double *rot, const double *vec);

double vector_dot(q_t n, const double *v1, const double *v2);

double *orthogonal_rotate(q_t n, const double *m1, const double *m2);

double *gen_rot(gsl_rng *r, q_t n);

