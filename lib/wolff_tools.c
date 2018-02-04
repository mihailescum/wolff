
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

