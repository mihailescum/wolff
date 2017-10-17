
#include "queue.h"
#include "wolff.h"

int32_t spin_to_sign(bool spin) {
  if (spin) {
    return -1;
  } else {
    return 1;
  }
}

int32_t sign(double x) {
  return x > 0 ? 1 : -1;
}

graph_t *graph_add_ext(const graph_t *g) {
  graph_t *h = (graph_t *)calloc(1, sizeof(graph_t));

  h->nv = g->nv + 1;
  h->ne = g->ne + g->nv;

  h->ev = (uint32_t *)malloc(2 * h->ne * sizeof(uint32_t));
  h->vei = (uint32_t *)malloc((h->nv + 1) * sizeof(uint32_t));
  h->ve = (uint32_t *)malloc(2 * h->ne * sizeof(uint32_t));
  h->vx = (double *)malloc(2 * h->nv * sizeof(double));
  h->bq = (bool *)malloc(h->nv * sizeof(bool));

  memcpy(h->ev, g->ev, 2 * g->ne * sizeof(uint32_t));
  memcpy(h->vx, g->vx, 2 * g->nv * sizeof(double));
  memcpy(h->bq, g->bq, g->nv * sizeof(bool));
  h->vx[2 * g->nv] = -1;
  h->vx[2 * g->nv + 1] = -0.5;
  h->bq[g->nv] = false;

  for (uint32_t i = 0; i < g->nv; i++) {
    h->ev[2 * g->ne + 2 * i] = i;
    h->ev[2 * g->ne + 2 * i + 1] = g->nv;
  }

  for (uint32_t i = 0; i < g->nv; i++) {
    h->vei[i] = g->vei[i] + i;

    for (uint32_t j = 0; j < g->vei[i + 1] - g->vei[i]; j++) {
      h->ve[h->vei[i] + j] = g->ve[g->vei[i] + j];
    }

    h->ve[h->vei[i] + g->vei[i + 1] - g->vei[i]] = g->ne + i;
  }

  h->vei[g->nv] = g->vei[g->nv] + g->nv;
  h->vei[g->nv + 1] = h->vei[g->nv] + g->nv;

  for (uint32_t i = 0; i < g->nv; i++) {
    h->ve[h->vei[g->nv] + i] = g->ne + i;
  }

  return h;
}

uint32_t get_neighbor(const graph_t *g, uint32_t v, uint32_t i) {
  uint32_t e, v1, v2;

  e = g->ve[g->vei[v] + i]; // select the ith bond connected to site
  v1 = g->ev[2 * e];
  v2 = g->ev[2 * e + 1];

  return v == v1 ? v2 : v1; // distinguish neighboring site from site itself
}

cluster_t *flip_cluster(const graph_t *g, const double *ps, bool *x, bool stop_on_ghost,
                        gsl_rng *r) {
  uint32_t v0;
  int32_t n_h_bonds, n_bonds;
  bool x0;
  cluster_t *c;

  v0 = gsl_rng_uniform_int(r, g->nv); // pick a random vertex
  x0 = x[v0];                         // record its orientation

  ll_t *stack = NULL;     // create a new stack
  stack_push(&stack, v0); // push the initial vertex to the stack

  // initiate the data structure for returning flip information
  c = (cluster_t *)calloc(1, sizeof(cluster_t));

  while (stack != NULL) {
    uint32_t v;
    uint32_t nn;

    v = stack_pop(&stack);
    nn = g->vei[v + 1] - g->vei[v];

    if (x[v] == x0 && !(c->hit_ghost && stop_on_ghost)) { // if the vertex hasn't already been flipped
      x[v] = !x[v];   // flip the vertex

      for (uint32_t i = 0; i < nn; i++) {
        bool is_ext;
        uint32_t vn;
        int32_t *bond_counter;
        double prob;

        vn = get_neighbor(g, v, i);

        is_ext = (v == g->nv - 1 || vn == g->nv - 1); // our edge contained the "ghost" spin if either of its vertices was the last in the graph

        bond_counter = is_ext ? &(c->dHb) : &(c->dJb);
        prob = is_ext ? ps[1] : ps[0];

        if (x[vn] ==
            x0) { // if the neighboring site matches the flipping cluster...
          (*bond_counter)++;

          if (gsl_rng_uniform(r) < prob) { // and with probability ps[e]...
            if (is_ext && stop_on_ghost) {
              c->hit_ghost = true;
            }
            stack_push(&stack, vn); // push the neighboring vertex to the stack
          }
        } else {
          (*bond_counter)--;
        }
      }

      if (v != g->nv - 1) { // count the number of non-external sites that flip
        c->nv++;
      }
    }
  }

  return c;
}

uint32_t wolff_step(double T, double H, ising_state_t *s, sim_t sim, gsl_rng *r,
                    double *ps) {
  uint32_t n_flips;
  bool no_r, no_ps;
  no_r = false;
  no_ps = false;

  if (r == NULL) {
    no_r = true;
    r = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(r, jst_rand_seed());
  }

  if (ps == NULL) {
    no_ps = true;
    ps = (double *)malloc(2 * sizeof(double));
    ps[0] = 1 - exp(-2 / T);
    ps[1] = 1 - exp(-2 * fabs(H) / T);
  }

  switch (sim) {
    case METROPOLIS: {
      uint32_t v0 = gsl_rng_uniform_int(r, s->g->nv); // pick a random vertex
      uint32_t nn = s->g->vei[v0 + 1] - s->g->vei[v0];

      double dE = 0;
      for (uint32_t i = 0; i < nn; i++) {
        bool is_ext;
        uint32_t vn;
        int32_t *bond_counter;
        double prob;

        vn = get_neighbor(s->g, v0, i);

        is_ext = (v0 == s->g->nv - 1 || vn == s->g->nv - 1); // our edge contained the "ghost" spin if either of its vertices was the last in the graph

        if (is_ext) {
          dE += 2 * H * spin_to_sign(s->spins[vn]) * spin_to_sign(s->spins[v0]);
        } else {
          dE += 2 * spin_to_sign(s->spins[vn]) * spin_to_sign(s->spins[v0]);
        }
      }

      double p = exp(-dE / T);
      if (gsl_rng_uniform(r) < p) {
        s->M += - sign(H) * 2 * spin_to_sign(s->spins[v0]) * spin_to_sign(s->spins[s->g->nv - 1]);
        s->H += dE;
        s->spins[v0] = !s->spins[v0];
      }

      n_flips = 1;
                     }
                     break;

    case WOLFF: {
      cluster_t *c = flip_cluster(s->g, ps, s->spins, false, r);
      s->M += - sign(H) * 2 * c->dHb;
      s->H += 2 * (c->dJb + sign (H) * H * c->dHb);
      n_flips = c->nv;

      free(c);
                }
                break;
    case WOLFF_GHOST: {
      bool *spins_bak;
      spins_bak = (bool *)malloc(s->g->nv * sizeof(bool));
      memcpy(spins_bak, s->spins, s->g->nv * sizeof(bool));

      cluster_t *c = flip_cluster(s->g, ps, s->spins, true, r);

      if (c->hit_ghost) {
        memcpy(s->spins, spins_bak, s->g->nv * sizeof(bool));
      } else {
        s->M += - sign(H) * 2 * c->dHb;
        s->H += 2 * (c->dJb + sign (H) * H * c->dHb);
      }
      free(spins_bak);

      n_flips = c->nv;

      free(c);
                      }
                break;
  }

  if (no_ps) {
    free(ps);
  }

  if (no_r) {
    gsl_rng_free(r);
  }

  return n_flips;
}

double add_to_avg(double mx, double x, uint64_t n) {
  return mx * (n / (n + 1.)) + x * 1. / (n + 1.);
}

void update_meas(meas_t *m, double x) {
  uint64_t n = m->n;

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
