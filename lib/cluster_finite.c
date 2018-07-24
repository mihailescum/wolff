
#include "cluster_finite.h"

v_t flip_cluster_finite(state_finite_t *s, v_t v0, R_t rot_ind, gsl_rng *r) {
  q_t *rot = s->transformations + s->q * s->involutions[rot_ind];
  q_t *R_inv = symmetric_invert(s->q, s->R);
  v_t nv = 0;

  ll_t *stack = NULL;     // create a new stack
  stack_push(&stack, v0); // push the initial vertex to the stack

  bool *marks = (bool *)calloc(s->g->nv, sizeof(bool));

  while (stack != NULL) {
    v_t v = stack_pop(&stack);

    if (!marks[v]) {
      q_t s_old, s_new;
      q_t *R_new, *R_inv_new; 
      bool external_flipped;

      marks[v] = true;

      if (v == s->g->nv - 1) {
        R_new = symmetric_compose(s->q, rot, s->R);
        R_inv_new = symmetric_invert(s->q, R_new);
        external_flipped = true;
      } else {
        s_old = s->spins[v];
        s_new = symmetric_act(rot, s_old);
        external_flipped = false;
      }

      v_t nn = s->g->v_i[v + 1] - s->g->v_i[v];

      for (v_t i = 0; i < nn; i++) {
        q_t sn, non_ghost;
        double prob;
        bool external_neighbor = false;

        v_t vn = s->g->v_adj[s->g->v_i[v] + i];

        if (vn == s->g->nv - 1) {
          external_neighbor = true;
          non_ghost = v;
        } else {
          sn = s->spins[vn];
          non_ghost = vn;
        }

        if (external_flipped || external_neighbor) {
          q_t rot_s_old, rot_s_new;

          if (external_neighbor) {
            rot_s_old = symmetric_act(R_inv, s_old);
            rot_s_new = symmetric_act(R_inv, s_new);
          } else {
            rot_s_old = symmetric_act(R_inv, sn);
            rot_s_new = symmetric_act(R_inv_new, sn);
          }

          prob = s->H_probs[rot_s_new * s->q + rot_s_old];

          s->M[rot_s_old]--;
          s->M[rot_s_new]++;

          for (D_t i = 0; i < s->D; i++) {
            L_t x = (non_ghost / (v_t)pow(s->L, s->D - i - 1)) % s->L;

            s->ReF[s->D * i + rot_s_old] -= s->precomputed_cos[i];
            s->ReF[s->D * i + rot_s_new] += s->precomputed_cos[i];

            s->ImF[s->D * i + rot_s_old] -= s->precomputed_sin[i];
            s->ImF[s->D * i + rot_s_new] += s->precomputed_sin[i];
          }

        } else {
          q_t diff_old = s->bond_with_zero_type[s->transformations[s->q * s->transform_site_to_zero[sn] + s_old]];
          q_t diff_new = s->bond_with_zero_type[s->transformations[s->q * s->transform_site_to_zero[sn] + s_new]];

          prob = s->J_probs[diff_new * s->n_bond_types + diff_old];

          s->B[diff_old]--;
          s->B[diff_new]++;
        }

        if (gsl_rng_uniform(r) < prob) { // and with probability ps[e]...
          stack_push(&stack, vn); // push the neighboring vertex to the stack
        }
      }

      if (external_flipped) {
        free(s->R);
        free(R_inv);
        s->R = R_new;
        R_inv = R_inv_new;
      } else {
        s->spins[v] = s_new;
      }

      if (v != s->g->nv - 1) { // count the number of non-external sites that flip
        nv++;
      }
    }
  }

  free(marks);
  free(R_inv);

  return nv;
}

