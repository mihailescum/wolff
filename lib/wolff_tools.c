
#include "wolff.h"

graph_t *graph_add_ext(const graph_t *g) {
  graph_t *h = (graph_t *)calloc(1, sizeof(graph_t));

  h->nv = g->nv + 1;
  h->ne = g->ne + g->nv;

  h->ev = (uint32_t *)malloc(2 * h->ne * sizeof(uint32_t));
  h->vei = (uint32_t *)malloc((h->nv + 1) * sizeof(uint32_t));
  h->ve = (uint32_t *) malloc(2 * h->ne * sizeof(uint32_t));
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

double get_hamiltonian(graph_t *g, double *coupling, bool *x) {
  double hamiltonian = 0;

  for (uint32_t i = 0; i < g->ne; i++) {
    uint32_t v1, v2;

    v1 = g->ev[2 * i];
    v2 = g->ev[2 * i + 1];

    if (x[v1] == x[v2]) {
      hamiltonian -= coupling[i];
    } else {
      hamiltonian += coupling[i];
    }
  }

  return hamiltonian;
}

cluster_t *flip_cluster(const graph_t *g, const double *ps, double H, bool *x, gsl_rng *r) {
  uint32_t v0;
  bool x0;
  cluster_t *c;
  
  v0 = gsl_rng_uniform_int(r, g->nv); // pick a random vertex
  x0 = x[v0]; // record its orientation

  ll_t *stack = NULL; // create a new stack
  stack_push(&stack, v0); // push the initial vertex to the stack

  // initiate the data structure for returning flip information
  c = (cluster_t *)calloc(1, sizeof(cluster_t));
  c->nv = 0;
  c->dH = 0;

  while (stack != NULL) {
    uint32_t v;
    uint16_t nn;

    v = stack_pop(&stack);
    nn = g->vei[v + 1] - g->vei[v];

    if (x[v] == x0) { // if the vertex hasn't already been flipped
      x[v] = !x[v]; // flip the vertex

      for (uint16_t i = 0; i < nn; i++) {
        uint32_t e, v1, v2, vn;

        e = g->ve[g->vei[v] + i]; // select the ith bond connected to site
        v1 = g->ev[2 * e];
        v2 = g->ev[2 * e + 1];

        vn = v == v1 ? v2 : v1; // distinguish neighboring site from site itself

        if (x[vn] == x0) { // if the neighboring site matches the flipping cluster...
          if (v1 == g->nv - 1 || v2 == g->nv - 1) {
            c->dH += H;
          } else {
            c->dH += 1;
          }

          if (gsl_rng_uniform(r) < ps[e]) { // and with probability ps[e]...
            stack_push(&stack, vn); // push the neighboring vertex to the stack
          }
        } else {
          if (v1 == g->nv - 1 || v2 == g->nv - 1) {
            c->dH -= H;
          } else {
            c->dH -= 1;
          }
        }
      }

      if (v != g->nv - 1) {
        c->nv++;
      }
    }
  }

  if (x0) {
    c->nv = -c->nv;
  }

  return c;
}

double hh(double th) {
  return (th - pow(th, 3) / 1.16951) * (1 - 0.222389 * pow(th, 2) - 0.043547 * pow(th, 4) - 0.014809 * pow(th, 6) - 0.007168 * pow(th, 8));
}

double *get_bond_probs(double T, double H, ising_state_t *s) {
  double p = 1 - exp(-2 / T);
  double q = 1 - exp(-2 * H / T);

  double *ps = (double *)malloc(s->g->ne * sizeof(double));

  for (uint32_t i = 0; i < s->g->ne; i++) {
    uint32_t v1, v2;
    v1 = s->g->ev[2 * i];
    v2 = s->g->ev[2 * i + 1];
    if (v1 == s->g->nv - 1 || v2 == s->g->nv - 1) {
      ps[i] = q;
    } else {
      ps[i] = p;
    }
  }

  return ps;
}

int32_t wolff_step(double T, double H, ising_state_t *s, gsl_rng *r, double *ps) {
  if (r == NULL) {
    r = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(r, jst_rand_seed());
  }

  if (ps == NULL) {
    ps = get_bond_probs(T, H, s);
  }

  cluster_t *c = flip_cluster(s->g, ps, H, s->spins, r);

  s->M += 2 * c->nv;
  s->H += 2 * c->dH;

  int32_t n_flipped = c->nv;

  free(c);

  return n_flipped;
}

