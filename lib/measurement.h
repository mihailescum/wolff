
#include <math.h>
#include <stdlib.h>

#include "types.h"
#include "stack.h"

typedef struct {
  uint64_t n;
  double x;
  double x2;
  double m2;
  double m4;
} meas_t;

typedef struct {
  uint64_t n;
  uint64_t W;
  double *OO;
  dll_t *Op;
  double O;
  double O2;
} autocorr_t;

void meas_update(meas_t *m, double x);

double meas_dx(const meas_t *m);

double meas_c(const meas_t *m);

double meas_dc(const meas_t *m);

void update_autocorr(autocorr_t *OO, double O);

double rho(const autocorr_t *o, uint64_t i);

