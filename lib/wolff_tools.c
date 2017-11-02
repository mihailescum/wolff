
#include "queue.h"
#include "wolff.h"

int8_t spin_to_sign(bool spin) {
  /* takes a spin (represented by a bool) and converts it to an integer (plus
   * or minus one). our convention takes true to be negative spin and false to
   * be positive spin
   */
  if (spin) {
    return -1;
  } else {
    return 1;
  }
}

int8_t sign(double x) {
  // the sign function. zero returns positive sign.
  return x >= 0 ? 1 : -1;
}

graph_t *graph_add_ext(const graph_t *G) {
  /* takes a graph object G and returns tG, the same graph with an extra vertex
   * that is connected to every other vertex.
   */
  graph_t *tG = (graph_t *)calloc(1, sizeof(graph_t));

  tG->nv = G->nv + 1;
  tG->ne = G->ne + G->nv;

  tG->ev = (uint32_t *)malloc(2 * tG->ne * sizeof(uint32_t));
  tG->vei = (uint32_t *)malloc((tG->nv + 1) * sizeof(uint32_t));
  tG->ve = (uint32_t *)malloc(2 * tG->ne * sizeof(uint32_t));
  tG->vx = (double *)malloc(2 * tG->nv * sizeof(double));
  tG->bq = (bool *)malloc(tG->nv * sizeof(bool));

  memcpy(tG->ev, G->ev, 2 * G->ne * sizeof(uint32_t));
  memcpy(tG->vx, G->vx, 2 * G->nv * sizeof(double));
  memcpy(tG->bq, G->bq, G->nv * sizeof(bool));

  tG->vx[2 * G->nv] = -1;
  tG->vx[2 * G->nv + 1] = -0.5;
  tG->bq[G->nv] = false;

  for (uint32_t i = 0; i < G->nv; i++) {
    tG->ev[2 * G->ne + 2 * i] = i;
    tG->ev[2 * G->ne + 2 * i + 1] = G->nv;
  }

  for (uint32_t i = 0; i < G->nv; i++) {
    tG->vei[i] = G->vei[i] + i;

    for (uint32_t j = 0; j < G->vei[i + 1] - G->vei[i]; j++) {
      tG->ve[tG->vei[i] + j] = G->ve[G->vei[i] + j];
    }

    tG->ve[tG->vei[i] + G->vei[i + 1] - G->vei[i]] = G->ne + i;
  }

  tG->vei[G->nv] = G->vei[G->nv] + G->nv;
  tG->vei[G->nv + 1] = tG->vei[G->nv] + G->nv;

  for (uint32_t i = 0; i < G->nv; i++) {
    tG->ve[tG->vei[G->nv] + i] = G->ne + i;
  }

  return tG;
}

uint32_t get_neighbor(const graph_t *g, uint32_t v, uint32_t i) {
  // returns the index of the ith neighbor of vertex v on graph g.
  assert(i < g->vei[v + 1] - g->vei[v]); // don't request a neighbor over the total number of neighbors!
  
  uint32_t e, v1, v2;

  e = g->ve[g->vei[v] + i]; // select the ith bond connected to site
  v1 = g->ev[2 * e];
  v2 = g->ev[2 * e + 1];

  return v == v1 ? v2 : v1; // distinguish neighboring site from site itself
}

cluster_t *flip_cluster(const graph_t *g, const double *ps, bool *s, bool stop_on_ghost,
                        gsl_rng *r) {
  /* flips a wolff cluster on the graph g with initial state s. s is always
   * changed by the routine. the probability of adding a normal bond to the
   * cluster is passed by ps[0], and the probability of adding a bond with the
   * external field spin to the cluster is passed by ps[1]. if stop_on_ghost is
   * true, adding the external field spin to the cluster stops execution of the
   * routine. flip_cluster returns an object of type cluster_t, which encodes
   * information about the flipped cluster.
   */
  uint32_t v0;
  int32_t n_h_bonds, n_bonds;
  bool s0;
  cluster_t *c;

  v0 = gsl_rng_uniform_int(r, g->nv); // pick a random vertex
  s0 = s[v0];                         // record its orientation

  ll_t *stack = NULL;     // create a new stack
  stack_push(&stack, v0); // push the initial vertex to the stack

  // initiate the data structure for returning flip information
  c = (cluster_t *)calloc(1, sizeof(cluster_t));

  while (stack != NULL) {
    uint32_t v;
    uint32_t nn;

    v = stack_pop(&stack);
    nn = g->vei[v + 1] - g->vei[v];

    if (s[v] == s0 && !(c->hit_ghost && stop_on_ghost)) { // if the vertex hasn't already been flipped
      s[v] = !s[v];   // flip the vertex
      if (stop_on_ghost) {
        stack_push(&(c->spins), v);
      }

      for (uint32_t i = 0; i < nn; i++) {
        bool is_ext;
        uint32_t vn;
        int32_t *bond_counter;
        double prob;

        vn = get_neighbor(g, v, i);

        is_ext = (v == g->nv - 1 || vn == g->nv - 1); // our edge contained the "ghost" spin if either of its vertices was the last in the graph

        bond_counter = is_ext ? &(c->dHb) : &(c->dJb);
        prob = is_ext ? ps[1] : ps[0];

        if (s[vn] ==
            s0) { // if the neighboring site matches the flipping cluster...
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

  if (ps == NULL) { // computing exponentials is relatively expensive, so calling functions are given the option to supply values calculated once in the entire runtime
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
      } break;

    case WOLFF: {
      cluster_t *c = flip_cluster(s->g, ps, s->spins, false, r);
      s->M += - sign(H) * 2 * c->dHb;
      s->H += 2 * (c->dJb + sign (H) * H * c->dHb);
      n_flips = c->nv;

      free(c);
      } break;
    case WOLFF_GHOST: {
      cluster_t *c = flip_cluster(s->g, ps, s->spins, true, r);

      if (c->hit_ghost) {
        while (c->spins != NULL) {
          uint32_t v = stack_pop(&(c->spins));
          s->spins[v] = !s->spins[v];
          // if we hit the external spin, undo the cluster flip
        }
      } else {
        while (c->spins != NULL) {
          stack_pop(&(c->spins));
          // we have to clear the memory on the stack anyway...
        }
        s->M += - sign(H) * 2 * c->dHb;
        s->H += 2 * (c->dJb + sign (H) * H * c->dHb);
      }

      n_flips = c->nv;

      free(c);
      } break;
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

void update_autocorr(autocorr_t *OO, double O) {
  OO->O = add_to_avg(OO->O, O, OO->n);
  OO->O2 = add_to_avg(OO->O2, pow(O, 2), OO->n);

  dll_t *Otmp = OO->Op;
  dll_t *Osave;
  uint64_t t = 0;

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

double rho(autocorr_t *o, uint64_t i) {
  return (o->OO[i] - pow(o->O, 2)) / (o->O2 - pow(o->O, 2));
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
