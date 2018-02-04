
#include "measurement.h"

double add_to_avg(double mx, double x, count_t n) {
  return mx * (n / (n + 1.0)) + x * 1.0 / (n + 1.0);
}

void update_meas(meas_t *m, double x) {
  count_t n = m->n;

  m->x = add_to_avg(m->x, x, n);
  m->x2 = add_to_avg(m->x2, pow(x, 2), n);

  m->m2 = add_to_avg(m->m2, pow(x - m->x, 2), n);
  m->m4 = add_to_avg(m->m4, pow(x - m->x, 4), n);

  if (n > 1) {
    double s2 = n / (n - 1.) * (m->x2 - pow(m->x, 2));
    m->dx = sqrt(s2 / n);
    m->c = s2;
    m->dc = sqrt((m->m4 - (n - 3.)/(n - 1.) * pow(m->m2, 2)) / n);
  }

  (m->n)++;
}

void update_autocorr(autocorr_t *OO, double O) {
  OO->O = add_to_avg(OO->O, O, OO->n);
  OO->O2 = add_to_avg(OO->O2, pow(O, 2), OO->n);

  dll_t *Otmp = OO->Op;
  dll_t *Osave;
  count_t t = 0;

  while (Otmp != NULL) {
    OO->OO[t] = add_to_avg(OO->OO[t], O * (Otmp->x), OO->n - t - 1);
    t++;
    if (t == OO->W - 1) {
      Osave = Otmp;
    }
    Otmp = Otmp->next;
  }

  if (t == OO->W) {
    if (OO->W == 1) {
      free(OO->Op);
      OO->Op = NULL;
    } else {
      free(Osave->next);
      Osave->next = NULL;
    }
  }

  stack_push_d(&(OO->Op), O);

  OO->n++;
}

double rho(autocorr_t *o, count_t i) {
  return (o->OO[i] - pow(o->O, 2)) / (o->O2 - pow(o->O, 2));
}

