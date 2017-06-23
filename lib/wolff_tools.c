
#include "queue.h"
#include "wolff.h"

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

cluster_t *flip_cluster(const graph_t *g, const double *ps, bool *x,
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

    if (x[v] == x0) { // if the vertex hasn't already been flipped
      x[v] = !x[v];   // flip the vertex

      for (uint32_t i = 0; i < nn; i++) {
        bool is_ext;
        uint32_t e, v1, v2, vn;
        int32_t *bond_counter;
        double prob;

        e = g->ve[g->vei[v] + i]; // select the ith bond connected to site
        v1 = g->ev[2 * e];
        v2 = g->ev[2 * e + 1];

        vn = v == v1 ? v2 : v1; // distinguish neighboring site from site itself

        is_ext = (v1 == g->nv - 1 || v2 == g->nv - 1);

        bond_counter = is_ext ? &(c->dHb) : &(c->dJb);
        prob = is_ext ? ps[1] : ps[0];

        if (x[vn] ==
            x0) { // if the neighboring site matches the flipping cluster...
          (*bond_counter)++;

          if (gsl_rng_uniform(r) < prob) { // and with probability ps[e]...
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

uint32_t wolff_step(double T, double H, ising_state_t *s, gsl_rng *r,
                    double *ps) {
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

  cluster_t *c = flip_cluster(s->g, ps, s->spins, r);

  s->M += - sign(H) * 2 * c->dHb;
  s->H += 2 * (c->dJb + sign (H) * H * c->dHb);

  uint32_t n_flips = c->nv;

  free(c);

  if (no_ps) {
    free(ps);
  }

  if (no_r) {
    gsl_rng_free(r);
  }

  return n_flips;
}
