
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

R_t *transformation_bringing_to_zero(q_t q, R_t n_transformations, q_t *transformations) {
  R_t *destination = (R_t *)malloc(q * sizeof(R_t));

  for (q_t i = 0; i < q; i++) {
    for (R_t j = 0; j < n_transformations; j++) {
      if (transformations[q * j + i] == 0) {
        destination[i] = j;
      }
    }
  }

  return destination;
}

R_t find_involutions(R_t *destination, q_t q, R_t n_transformations, q_t *transformations) {
  R_t n_involutions = 0;

  for (R_t i = 1; i < n_transformations; i++) {
    bool is_involution = true;
    for (q_t j = 0; j < q; j++) {
      if (j != transformations[q * i + transformations[q * i + j]]) {
        is_involution = false;
        break;
      }
    }
    if (is_involution) {
      destination[n_involutions] = i;
      n_involutions++;
    }
  }

  return n_involutions;
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

  s->n_transformations = 2;
  s->transformations = (q_t *)malloc(2 * 2 * sizeof(q_t));
  s->transformations[0] = 0;
  s->transformations[1] = 1;
  s->transformations[2] = 1;
  s->transformations[3] = 0;

  s->n_involutions = 1;
  s->involutions = (R_t *)malloc(1 * sizeof(R_t));
  s->involutions[0] = 1;

  s->transform_site_to_zero = (R_t *)malloc(2 * sizeof(R_t));
  s->transform_site_to_zero[0] = 0;
  s->transform_site_to_zero[1] = 1;

  s->n_bond_types = 2;
  s->bond_with_zero_type = (q_t *)malloc(2 * sizeof(q_t));
  s->bond_with_zero_type[0] = 0;
  s->bond_with_zero_type[1] = 1;

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

  s->M = (v_t *)calloc(2, sizeof(v_t));
  s->M[0] = s->nv; // everyone starts in state 0, remember?
  s->B = (v_t *)calloc(2, sizeof(v_t));
  s->B[0] = s->ne;

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
  s->n_transformations = factorial(q);
  s->transformations = symmetric_gen_transformations(q);
  s->involutions = (R_t *)malloc(s->n_transformations * sizeof(R_t));
  s->n_involutions = find_involutions(s->involutions, q, s->n_transformations, s->transformations);

  s->transform_site_to_zero = transformation_bringing_to_zero(q, s->n_transformations, s->transformations);

  s->n_bond_types = 2;

  s->bond_with_zero_type = (q_t *)malloc(q * sizeof(q_t));

  s->bond_with_zero_type[0] = 0;

  for (q_t i = 1; i < q; i++) {
    s->bond_with_zero_type[i] = 1;
  }

  s->T = T;
  s->J = (double *)calloc(2, sizeof(double)); 
  s->J[0] = 1.0;
  s->J[1] = 0.0;

  s->H = (double *)malloc(q * sizeof(double)); 
  for (q_t i = 0; i < q; i++) {
    s->H[i] = H[i];
  }

  s->J_probs = Jprobs_from_J(s->n_bond_types, T, s->J);
  s->H_probs = Jprobs_from_J(q, T, s->H);

  s->spins = (q_t *)calloc(s->nv, sizeof(q_t));
  s->R = initialize_R(q);

  s->M = (v_t *)calloc(q, sizeof(v_t));
  s->M[0] = s->nv; // everyone starts in state 0, remember?
  s->B = (v_t *)calloc(s->n_bond_types, sizeof(v_t));
  s->B[0] = s->ne; // everyone starts in state 0, remember?

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

  s->n_transformations = 2 * q;
  s->transformations = dihedral_gen_transformations(q);
  s->n_involutions = q;
  s->involutions = dihedral_gen_involutions(q);

  s->transform_site_to_zero = transformation_bringing_to_zero(q, s->n_transformations, s->transformations);
  s->bond_with_zero_type = malloc(q * sizeof(q_t));

  s->n_bond_types = q / 2 + 1;

  for (q_t i = 0; i < q / 2 + 1; i++) {
    s->bond_with_zero_type[i] = i;
  }

  for (q_t i = 1; i < (q + 1) / 2; i++) { 
    s->bond_with_zero_type[q - i] = i;
  }


  s->T = T;
  s->J = (double *)malloc(s->n_bond_types * sizeof(double)); 

  for (q_t i = 0; i < s->n_bond_types; i++) {
    s->J[i] = cos(2 * M_PI * i / ((double)q));
  }


  s->H = (double *)malloc(q * sizeof(double)); 
  for (q_t i = 0; i < q; i++) {
    s->H[i] = H[i];
  }

  s->J_probs = Jprobs_from_J(s->n_bond_types, T, s->J);
  s->H_probs = Jprobs_from_J(q, T, s->H);

  s->spins = (q_t *)calloc(s->nv, sizeof(q_t));
  s->R = initialize_R(q);

  s->M = (v_t *)calloc(q, sizeof(v_t));
  s->M[0] = s->nv; // everyone starts in state 0, remember?
  s->B = (v_t *)calloc(s->n_bond_types, sizeof(v_t));
  s->B[0] = s->ne; // everyone starts in state 0, remember?

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

  s->n_transformations = 2 * q;
  s->transformations = dihedral_gen_transformations(q);
  s->n_involutions = q;
  s->involutions = dihedral_gen_involutions(q);

  s->transform_site_to_zero = transformation_bringing_to_zero(q, s->n_transformations, s->transformations);
  s->bond_with_zero_type = malloc(q * sizeof(q_t));

  s->n_bond_types = q / 2 + 1;

  for (q_t i = 0; i < q / 2 + 1; i++) {
    s->bond_with_zero_type[i] = i;
  }

  for (q_t i = 1; i < (q + 1) / 2; i++) { 
    s->bond_with_zero_type[(int)q - (int)i] = i;
  }

  s->T = T;
  s->J = (double *)malloc(s->n_bond_types * sizeof(double)); 

  for (q_t i = 0; i < s->n_bond_types; i++) {
    s->J[i] = -pow(i, 2);
  }

  s->H = (double *)malloc(q * sizeof(double)); 
  for (q_t i = 0; i < q; i++) {
    s->H[i] = H[i];
  }

  s->J_probs = Jprobs_from_J(s->n_bond_types, T, s->J);
  s->H_probs = Jprobs_from_J(q, T, s->H);

  s->spins = (q_t *)calloc(s->nv, sizeof(q_t));
  s->R = initialize_R(q);

  s->M = (v_t *)calloc(q, sizeof(v_t));
  s->M[0] = s->nv; // everyone starts in state 0, remember?
  s->B = (v_t *)calloc(s->n_bond_types, sizeof(v_t));
  s->B[0] = s->nv; // everyone starts in state 0, remember?

  return s;
}

double state_finite_energy(state_finite_t *s) {
  double E = 0;

  for (q_t i = 0; i < s->n_bond_types; i++) {
    E += s->J[i] * s->B[i];
  }
  for (q_t i = 0; i < s->q; i++) {
    E += s->H[i] * s->M[i];
  }

  return -E;
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
  free(s->B);
  free(s->transformations);
  free(s->involutions);
  free(s->transform_site_to_zero);
  free(s->bond_with_zero_type);
  free(s);
}

