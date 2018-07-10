
#include "convex.h"
#include "measurement.h"

meas_t *meas_initialize(count_t W) {
  meas_t *m = (meas_t *)calloc(1, sizeof(meas_t));
  m->W = W;
  m->xx = (double *)calloc(2 * W + 1, sizeof(double));

  return m;
}

double add_to_avg(double mx, double x, count_t n) {
  return mx * (n / (n + 1.0)) + x / (n + 1.0);
}

void meas_update(meas_t *m, double x) {
  count_t n = m->n;

  m->x = add_to_avg(m->x, x, n);
  m->x2 = add_to_avg(m->x2, pow(x, 2), n);
  m->x4 = add_to_avg(m->x4, pow(x, 4), n);

  m->m2 = add_to_avg(m->m2, pow(x - m->x, 2), n);
  m->m4 = add_to_avg(m->m4, pow(x - m->x, 4), n);

  dll_t *tmp_window = m->x_window;
  dll_t *pos_save;
  count_t t = 0;

  while (tmp_window != NULL) {
    m->xx[t] = add_to_avg(m->xx[t], x * (tmp_window->x), m->n - t - 1);
    t++;
    if (t == 2 * (m->W)) {
      pos_save = tmp_window;
    }
    tmp_window = tmp_window->next;
  }

  if (t == 2 * (m->W) + 1) {
    if (2 * (m->W) + 1 == 1) {
      free(m->x_window);
      m->x_window = NULL;
    } else {
      free(pos_save->next);
      pos_save->next = NULL;
    }
  }

  stack_push_d(&(m->x_window), x);

  (m->n)++;
}

double meas_dx(const meas_t *m) {
  return 2 * get_tau(m) * Cxx(m, 0) / m->n;
//  return sqrt(1. / (m->n - 1.) * (m->x2 - pow(m->x, 2)));
}

double meas_c(const meas_t *m) {
  return m->n / (m->n - 1.) * (m->x2 - pow(m->x, 2));
}

double meas_dc(const meas_t *m) {
  return sqrt((m->m4 - (m->n - 3.)/(m->n - 1.) * pow(m->m2, 2)) / m->n);
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

double rho(const autocorr_t *o, count_t i) {
  return (o->OO[i] - pow(o->O, 2)) / (o->O2 - pow(o->O, 2));
}

double Cxx(const meas_t *m, count_t t) {
  return m->xx[t] - pow(m->x, 2);
}

double rho_m(const meas_t *m, count_t t) {
  return Cxx(m, t) / Cxx(m, 0);
}

double get_tau(const meas_t *m) {
  double *Gammas = (double *)malloc((m->W + 1) * sizeof(double));

  Gammas[0] = 1 + rho_m(m, 0);
  for (uint64_t i = 0; i < m->W; i++) {
    Gammas[1 + i] = rho_m(m, 2 * i + 1) + rho_m(m, 2 * i + 2);
  }

  uint64_t n;
  for (n = 0; n < m->W + 1; n++) {
    if (Gammas[n] <= 0) {
      break;
    }
  }

  double *conv_Gamma = get_convex_minorant(n, Gammas);

  double tau = - 0.5;

  for (uint64_t i = 0; i < n + 1; i++) {
    tau += conv_Gamma[i];
  }

  free(Gammas);

  return tau;
}

void print_meas(const meas_t *m, const char *sym, FILE *outfile) {
  fprintf(outfile, "%s-><|n->%" PRIcount ",x->%.15f,x^2->%.15f,x^4->%.15f,xx->{", sym, m->n, m->x, m->x2, m->x4);
  for (count_t i = 0; i < 2 * (m->W) + 1; i++) {
    fprintf(outfile, "%.15f", m->xx[i]);
    if (i < 2 * (m->W)) {
      fprintf(outfile, ",");
    }
  }
  fprintf(outfile, "}|>");
}

void print_vec_meas(q_t q, const meas_t **m, const char *sym, FILE *outfile) {
  fprintf(outfile, "%s-><|n->{", sym);
  for (q_t i = 0; i < q; i++) {
    fprintf(outfile, "%" PRIcount, m[i]->n);
    if (i < q - 1) {
      fprintf(outfile, ",");
    }
  }
  fprintf(outfile, "},x->{");
  for (q_t i = 0; i < q; i++) {
    fprintf(outfile, "%.15f", m[i]->x);
    if (i < q - 1) {
      fprintf(outfile, ",");
    }
  }
  fprintf(outfile, "},x^2->{");
  for (q_t i = 0; i < q; i++) {
    fprintf(outfile, "%.15f", m[i]->x2);
    if (i < q - 1) {
      fprintf(outfile, ",");
    }
  }
  fprintf(outfile, "},x^4->{");
  for (q_t i = 0; i < q; i++) {
    fprintf(outfile, "%.15f", m[i]->x4);
    if (i < q - 1) {
      fprintf(outfile, ",");
    }
  }
  fprintf(outfile, "},xx->{");
  for (q_t i = 0; i < q; i++) {
    fprintf(outfile, "{");
    for (count_t j = 0; j < 2 * (m[i]->W) + 1; j++) {
      fprintf(outfile, "%.15f", m[i]->xx[j]);
      if (j < 2 * (m[i]->W)) {
        fprintf(outfile, ",");
      }
    }
    fprintf(outfile, "}");
    if (i < q - 1) {
      fprintf(outfile, ",");
    }
  }
  fprintf(outfile, "}|>");
}

void free_meas(meas_t *m) {
  free(m->xx);
  while (m->x_window != NULL) {
    stack_pop_d(&(m->x_window));
  }
  free(m);
}


