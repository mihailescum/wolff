
#include "initial_finite.h"

double *Jprobs_from_J(q_t q, double T, double *J) {
  double *J_probs = (double *)calloc(pow(q, 2), sizeof(double));

  for (q_t i = 0; i < q; i++) {
    for (q_t j = 0; j < q; j++) {
      J_probs[q * i + j] = 1.0 - exp((J[i] - J[j]) / T);
    }
  }

  return J_probs;
}

q_t *initialize_R(q_t q) {
  q_t *R = (q_t *)malloc(q * sizeof(q_t));

  for (q_t i = 0; i < q; i++) {
    R[i] = i;
  }

  return R;
}

state_finite_t *initial_finite_prepare_ising(D_t D, L_t L, double T, double *H) {
  state_finite_t *s = (state_finite_t *)calloc(1, sizeof(state_finite_t));

  s->D = D;
  s->L = L;

  {
    graph_t *g = graph_create_square(D, L);
    s->nv = g->nv;
    s->ne = g->ne;
    s->g = graph_add_ext(g);
    graph_free(g);
  }

  s->q = 2;
  s->n_transformations = 1;

  s->transformations = (q_t *)malloc(2 * sizeof(q_t));
  s->transformations[0] = 1;
  s->transformations[1] = 0;

  s->T = T;
  s->J = (double *)malloc(2 * sizeof(double)); 
  s->J[0] = 1.0;
  s->J[1] = -1.0;
  s->H = (double *)malloc(2 * sizeof(double)); 
  s->H[0] = H[0];
  s->H[1] = -H[0];

  s->J_probs = Jprobs_from_J(2, T, s->J);
  s->H_probs = Jprobs_from_J(2, T, s->H);

  s->spins = (q_t *)calloc(s->nv, sizeof(q_t));
  s->R = initialize_R(2);

  s->E = - ((double)s->ne) * s->J[0] - ((double)s->nv) * s->H[0];
  s->M = (v_t *)calloc(2, sizeof(v_t));
  s->M[0] = s->nv; // everyone starts in state 0, remember?

  return s;
}

state_finite_t *initial_finite_prepare_potts(D_t D, L_t L, q_t q, double T, double *H) {
  state_finite_t *s = (state_finite_t *)calloc(1, sizeof(state_finite_t));

  s->D = D;
  s->L = L;

  {
    graph_t *g = graph_create_square(D, L);
    s->nv = g->nv;
    s->ne = g->ne;
    s->g = graph_add_ext(g);
    graph_free(g);
  }

  s->q = q;
  s->n_transformations = q;
  s->transformations = dihedral_gen_transformations(q);

  s->T = T;
  s->J = (double *)calloc(q, sizeof(double)); 
  s->J[0] = 1.0;

  s->H = (double *)malloc(q * sizeof(double)); 
  for (q_t i = 0; i < q; i++) {
    s->H[i] = H[i];
  }

  s->J_probs = Jprobs_from_J(q, T, s->J);
  s->H_probs = Jprobs_from_J(q, T, s->H);

  s->spins = (q_t *)calloc(s->nv, sizeof(q_t));
  s->R = initialize_R(q);

  s->E = - ((double)s->ne) * s->J[0] - ((double)s->nv) * s->H[0];
  s->M = (v_t *)calloc(q, sizeof(v_t));
  s->M[0] = s->nv; // everyone starts in state 0, remember?

  return s;
}

state_finite_t *initial_finite_prepare_clock(D_t D, L_t L, q_t q, double T, double *H) {
  state_finite_t *s = (state_finite_t *)calloc(1, sizeof(state_finite_t));

  s->D = D;
  s->L = L;

  {
    graph_t *g = graph_create_square(D, L);
    s->nv = g->nv;
    s->ne = g->ne;
    s->g = graph_add_ext(g);
    graph_free(g);
  }

  s->q = q;
  s->n_transformations = q;
  s->transformations = dihedral_gen_transformations(q);

  s->T = T;
  s->J = (double *)malloc(q * sizeof(double)); 

  for (q_t i = 0; i < q; i++) {
    s->J[i] = cos(2 * M_PI * i / ((double)q));
  }


  s->H = (double *)malloc(q * sizeof(double)); 
  for (q_t i = 0; i < q; i++) {
    s->H[i] = H[i];
  }

  s->J_probs = Jprobs_from_J(q, T, s->J);
  s->H_probs = Jprobs_from_J(q, T, s->H);

  s->spins = (q_t *)calloc(s->nv, sizeof(q_t));
  s->R = initialize_R(q);

  s->E = - ((double)s->ne) * s->J[0] - ((double)s->nv) * s->H[0];
  s->M = (v_t *)calloc(q, sizeof(v_t));
  s->M[0] = s->nv; // everyone starts in state 0, remember?

  return s;
}


state_finite_t *initial_finite_prepare_dgm(D_t D, L_t L, q_t q, double T, double *H) {
  state_finite_t *s = (state_finite_t *)calloc(1, sizeof(state_finite_t));

  s->D = D;
  s->L = L;

  {
    graph_t *g = graph_create_square(D, L);
    s->nv = g->nv;
    s->ne = g->ne;
    s->g = graph_add_ext(g);
    graph_free(g);
  }

  s->q = q;
  s->n_transformations = q;
  s->transformations = dihedral_gen_transformations(q);

  s->T = T;
  s->J = (double *)malloc(q * sizeof(double)); 

  for (q_t i = 0; i < q / 2 + 1; i++) {
    s->J[i] = -pow(i, 2);
  }
  for (q_t i = 1; i < (q + 1) / 2; i++) {
    s->J[q - i] = -pow(i, 2);
  }

  s->H = (double *)malloc(q * sizeof(double)); 
  for (q_t i = 0; i < q; i++) {
    s->H[i] = H[i];
  }

  s->J_probs = Jprobs_from_J(q, T, s->J);
  s->H_probs = Jprobs_from_J(q, T, s->H);

  s->spins = (q_t *)calloc(s->nv, sizeof(q_t));
  s->R = initialize_R(q);

  s->E = - ((double)s->ne) * s->J[0] - ((double)s->nv) * s->H[0];
  s->M = (v_t *)calloc(q, sizeof(v_t));
  s->M[0] = s->nv; // everyone starts in state 0, remember?

  return s;
}

void state_finite_free(state_finite_t *s) {
  graph_free(s->g);
  free(s->J);
  free(s->H);
  free(s->J_probs);
  free(s->H_probs);
  free(s->spins);
  free(s->R);
  free(s->M);
  free(s->transformations);
  free(s);
}

