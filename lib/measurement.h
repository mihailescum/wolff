
#include <math.h>
#include <stdlib.h>

#include "types.h"
#include "stack.h"

typedef struct {
  uint64_t n;
  double x;
  double dx;
  double x2;
  double m2;
  double m4;
  double c;
  double dc;
} meas_t;

typedef struct {
  uint64_t n;
  uint64_t W;
  double *OO;
  dll_t *Op;
  double O;
  double O2;
} autocorr_t;

void update_meas(meas_t *m, double x);

void update_autocorr(autocorr_t *OO, double O);

double rho(autocorr_t *o, uint64_t i);

