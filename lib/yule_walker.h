
#pragma once

#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

#include "types.h"
#include "measurement.h"

double yule_walker(const autocorr_t *a);

