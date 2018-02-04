
#include "wolff.h"

v_t flip_cluster(ising_state_t *s, v_t v0, q_t step, gsl_rng *r) {
  q_t s0 = s->spins[v0];
  v_t nv = 0;

  ll_t *stack = NULL;     // create a new stack
  stack_push(&stack, v0); // push the initial vertex to the stack

  //node_t *T = NULL;
  bool *marks = (bool *)calloc(s->g->nv, sizeof(bool));

  while (stack != NULL) {
    v_t v = stack_pop(&stack);

//    if (!tree_contains(T, v)) { // if the vertex hasn't already been flipped
    if (!marks[v]) {
      q_t s_old = s->spins[v];
      q_t s_new = (s->spins[v] + step) % s->q;

      s->spins[v] = s_new;   // flip the vertex
      //tree_insert(&T, v);
      marks[v] = true;

      v_t nn = s->g->v_i[v + 1] - s->g->v_i[v];

      for (v_t i = 0; i < nn; i++) {
        v_t vn = s->g->v_adj[s->g->v_i[v] + i];
        q_t sn = s->spins[vn];
        double prob;

        bool is_ext = (v == s->g->nv - 1 || vn == s->g->nv - 1);

        q_t M_ind_0;
        q_t M_ind_1;

        if (is_ext) {
          if (vn == s->g->nv - 1) {
            M_ind_0 = (s_old + s->q - sn) % s->q;
            M_ind_1 = (s_new + s->q - sn) % s->q;
          } else {
            M_ind_0 = (sn + s->q - s_old) % s->q;
            M_ind_1 = (sn + s->q - s_new) % s->q;
          }
          prob = s->H_probs[M_ind_1 * s->q + M_ind_0];
          s->M[M_ind_0]--;
          s->M[M_ind_1]++;
          s->E += - s->H[M_ind_1] + s->H[M_ind_0];
        } else {
          M_ind_0 = (s_old + s->q - sn) % s->q;
          M_ind_1 = (s_new + s->q - sn) % s->q;
          prob = s->J_probs[M_ind_1 * s->q + M_ind_0];
          s->E += - s->J[M_ind_1] + s->J[M_ind_0];
        }

        if (gsl_rng_uniform(r) < prob) { // and with probability ps[e]...
          stack_push(&stack, vn); // push the neighboring vertex to the stack
        }
      }

      if (v != s->g->nv - 1) { // count the number of non-external sites that flip
        nv++;
      }
    }
  }

  //tree_freeNode(T);
  free(marks);

  return nv;
}

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
