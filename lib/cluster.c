
#include "cluster.h"

v_t flip_cluster(ising_state_t *s, v_t v0, q_t rot, gsl_rng *r) {
  v_t nv = 0;

  ll_t *stack = NULL;     // create a new stack
  stack_push(&stack, v0); // push the initial vertex to the stack

  bool *marks = (bool *)calloc(s->g->nv, sizeof(bool));

  while (stack != NULL) {
    v_t v = stack_pop(&stack);

    if (!marks[v]) {
      q_t s_old, s_new;
      dihedral_t *R_new; 
      bool external_flipped;

      marks[v] = true;

      if (v == s->g->nv - 1) {
        R_new = dihedral_compose(s->q, rot, s->R);
        external_flipped = true;
      } else {
        s_old = s->spins[v];
        s_new = dihedral_act(s->q, rot, s_old);
        external_flipped = false;
      }

      v_t nn = s->g->v_i[v + 1] - s->g->v_i[v];

      for (v_t i = 0; i < nn; i++) {
        q_t sn;
        double prob;
        bool external_neighbor = false;

        v_t vn = s->g->v_adj[s->g->v_i[v] + i];

        if (vn == s->g->nv - 1) {
          external_neighbor = true;
        } else {
          sn = s->spins[vn];
        }

        if (external_flipped || external_neighbor) {
          q_t rot_s_old, rot_s_new;

          if (external_neighbor) {
            rot_s_old = dihedral_inverse_act(s->q, s->R, s_old);
            rot_s_new = dihedral_inverse_act(s->q, s->R, s_new);
          } else {
            rot_s_old = dihedral_inverse_act(s->q, s->R, sn);
            rot_s_new = dihedral_inverse_act(s->q, R_new, sn);
          }

          prob = s->H_probs[rot_s_new * s->q + rot_s_old];

          s->M[rot_s_old]--;
          s->M[rot_s_new]++;

          s->E += - s->H[rot_s_new] + s->H[rot_s_old];
        } else {
          q_t diff_old = (s_old + s->q - sn) % s->q;
          q_t diff_new = (s_new + s->q - sn) % s->q;

          prob = s->J_probs[diff_new * s->q + diff_old];

          s->E += - s->J[diff_new] + s->J[diff_old];
        }

        if (gsl_rng_uniform(r) < prob) { // and with probability ps[e]...
          stack_push(&stack, vn); // push the neighboring vertex to the stack
        }
      }

      if (external_flipped) {
        free(s->R);
        s->R = R_new;
      } else {
        s->spins[v] = s_new;
      }

      if (v != s->g->nv - 1) { // count the number of non-external sites that flip
        nv++;
      }
    }
  }

  free(marks);

  return nv;
}

v_t flip_cluster_dgm(dgm_state_t *s, v_t v0, h_t rot, gsl_rng *r) {
  v_t nv = 0;

  ll_t *stack = NULL;     // create a new stack
  stack_push(&stack, v0); // push the initial vertex to the stack

  bool *marks = (bool *)calloc(s->g->nv, sizeof(bool));

  while (stack != NULL) {
    v_t v = stack_pop(&stack);

    if (!marks[v]) {
      h_t s_old, s_new;
      dihinf_t *R_new; 
      bool external_flipped;

      marks[v] = true;

      if (v == s->g->nv - 1) {
        R_new = dihinf_compose(rot, s->R);
        external_flipped = true;
      } else {
        s_old = s->spins[v];
        s_new = dihinf_act(rot, s_old);
        external_flipped = false;
      }

      v_t nn = s->g->v_i[v + 1] - s->g->v_i[v];

      for (v_t i = 0; i < nn; i++) {
        h_t sn;
        double prob;
        bool external_neighbor = false;

        v_t vn = s->g->v_adj[s->g->v_i[v] + i];

        if (vn == s->g->nv - 1) {
          external_neighbor = true;
        } else {
          sn = s->spins[vn];
        }

        if (external_flipped || external_neighbor) {
          h_t rot_s_old, rot_s_new;

          if (external_neighbor) {
            rot_s_old = dihinf_inverse_act(s->R, s_old);
            rot_s_new = dihinf_inverse_act(s->R, s_new);
          } else {
            rot_s_old = dihinf_inverse_act(s->R, sn);
            rot_s_new = dihinf_inverse_act(R_new, sn);
          }

          double dE = s->H(s->H_info, rot_s_old) - s->H(s->H_info, rot_s_new);
          prob = 1.0 - exp(-dE / s->T);

          s->M += rot_s_new - rot_s_old;
          s->E += dE;
        } else {
          double dE = (s->J)(s_old - sn) - (s->J)(s_new - sn);
          prob = 1.0 - exp(-dE / s->T);
          s->E += dE;
        }

        if (gsl_rng_uniform(r) < prob) { // and with probability ps[e]...
          stack_push(&stack, vn); // push the neighboring vertex to the stack
        }
      }

      if (external_flipped) {
        free(s->R);
        s->R = R_new;
      } else {
        s->spins[v] = s_new;
      }

      if (v != s->g->nv - 1) { // count the number of non-external sites that flip
        nv++;
      }
    }
  }

  free(marks);

  return nv;
}

v_t flip_cluster_vector(vector_state_t *s, v_t v0, double *rot, gsl_rng *r) {
  v_t nv = 0;

  ll_t *stack = NULL;     // create a new stack
  stack_push(&stack, v0); // push the initial vertex to the stack

  //node_t *T = NULL;
  bool *marks = (bool *)calloc(s->g->nv, sizeof(bool));

  while (stack != NULL) {
    v_t v = stack_pop(&stack);

//    if (!tree_contains(T, v)) { // if the vertex hasn't already been flipped
    if (!marks[v]) {
      double *s_old, *s_new, *R_tmp; 

      //tree_insert(&T, v);
      marks[v] = true;

      if (v == s->g->nv - 1) {
        R_tmp = orthogonal_rotate(s->n, rot, s->R);
      } else {
        s_old = &(s->spins[s->n * v]); // don't free me! I'm a pointer within array s->spins
        s_new = vector_rotate(s->n, rot, s_old); // free me! I'm a new vector
      }

      v_t nn = s->g->v_i[v + 1] - s->g->v_i[v];

      for (v_t i = 0; i < nn; i++) {
        v_t vn = s->g->v_adj[s->g->v_i[v] + i];
        double *sn;
        if (vn != s->g->nv - 1) {
          sn = &(s->spins[s->n * vn]);
        }
        double prob;

        bool is_ext = (v == s->g->nv - 1 || vn == s->g->nv - 1);

        if (is_ext) {
          double *rs_old, *rs_new;
          if (vn == s->g->nv - 1) {
            rs_old = vector_rotate_inverse(s->n, s->R, s_old);
            rs_new = vector_rotate_inverse(s->n, s->R, s_new);
          } else {
            rs_old = vector_rotate_inverse(s->n, s->R, sn);
            rs_new = vector_rotate_inverse(s->n, R_tmp, sn);
          }
          double dE = s->H(s->n, s->H_info, rs_old) - s->H(s->n, s->H_info, rs_new);
          prob = 1.0 - exp(-dE / s->T);
          vector_subtract(s->n, s->M, rs_old);
          vector_add(s->n, s->M, rs_new);
          s->E += dE;

          free(rs_old);
          free(rs_new);
        } else {
          double dE = (s->J)(vector_dot(s->n, sn, s_old)) - (s->J)(vector_dot(s->n, sn, s_new));
          prob = 1.0 - exp(-dE / s->T);
          //printf("(%g %g) (%g %g) (%g %g) %g\n", s_old[0], s_old[1], s_new[0], s_new[1], sn[0], sn[1], dE);
          //getchar();
          s->E += dE;
        }

        if (gsl_rng_uniform(r) < prob) { // and with probability ps[e]...
          stack_push(&stack, vn); // push the neighboring vertex to the stack
        }
      }

      if (v == s->g->nv - 1) {
        free(s->R);
        s->R = R_tmp;
      } else {
        vector_replace(s->n, s_old, s_new);
        free(s_new);
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

/*G
template <class R_t, class X_t>
v_t flip_cluster(state_t <R_t, X_t> *state, v_t v0, R_t *r, gsl_rng *rand) {
  v_t nv = 0;

  ll_t *stack = NULL;     // create a new stack
  stack_push(&stack, v0); // push the initial vertex to the stack

  bool *marks = (bool *)calloc(state->g->nv, sizeof(bool));

  while (stack != NULL) {
    v_t v = stack_pop(&stack);

    if (!marks[v]) {
      X_t *si_old, *si_new;
      R_t *s0_old, *s0_new;

      si_old = state->s[v];
      s0_old = state->s0;

      marks[v] = true;

      if (v == state->g->nv - 1) {
        s0_new = act <R_t, R_t> (r, s0_old);
      } else {
        si_new = act <R_t, X_t> (r, si_old);
      }

      v_t nn = state->g->v_i[v + 1] - state->g->v_i[v];

      for (v_t i = 0; i < nn; i++) {
        v_t vn = state->g->v_adj[state->g->v_i[v] + i];

        X_t *sj;

        if (vn != state->g->nv - 1) {
          sj = state->s[vn];
        }

        double prob;

        bool is_ext = (v == state->g->nv - 1 || vn == state->g->nv - 1);

        if (is_ext) {
          X_t *rs_old, *rs_new;
          if (vn == state->g->nv - 1) {
            rs_old = inverse_act <class R_t, class X_t> (s0_old, si_old);
            rs_new = inverse_act <class R_t, class X_t> (s0_old, si_new);
          } else {
            rs_old = inverse_act <class R_t, class X_t> (s0_old, sj);
            rs_new = inverse_act <class R_t, class X_t> (s0_new, sj);
          }
          double dE = state->B(rs_old) - state->B(rs_new);
          prob = 1.0 - exp(-dE / state->T);
          update_magnetization <X_t> (state->M, rs_old, rs_new);
          state->E += dE;

          free_X <X_t> (rs_old);
          free_X <X_t> (rs_new);
        } else {
          double dE = state->Z(si_old, sj) - state->Z(si_new, sj);
          prob = 1.0 - exp(-dE / state->T);
          state->E += dE;
        }

        if (gsl_rng_uniform(rand) < prob) { // and with probability...
          stack_push(&stack, vn); // push the neighboring vertex to the stack
        }
      }

      if (v == state->g->nv - 1) {
        free_R <R_t> (state->s0);
        state->s0 = s0_new;
      } else {
        free_X <X_t> (state->s[v]);
        state->s[v] = si_new;
      }

      if (v != state->g->nv - 1) { // count the number of non-external sites that flip
        nv++;
      }
    }
  }

  free(marks);

  return nv;
}
*/
