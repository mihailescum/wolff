
#pragma once

#include <functional>
#include <assert.h>
#include <fftw3.h>
#include <float.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <inttypes.h>
#include <cmath>
#include <stdbool.h>
#include <string.h>
#include <sys/types.h>

#include "types.h"
#include "rand.h"
#include "stack.h"
#include "convex.h"
#include "graph.h"
#include "tree.h"
#include "measurement.h"
#include "orthogonal.h"
#include "dihedral.h"
#include "dihinf.h"
#include "yule_walker.h"

template <class T>
void init(T*);

template <class T>
T scalar_multiple(v_t a, T b);

template <class R_t, class X_t>
X_t act(R_t a, X_t b);

template <class R_t, class X_t>
X_t act_inverse(R_t a, X_t b);

template <class T>
T copy(T a);

template <class T>
void free_spin(T a);

template <class T>
T add(T, T);

template <class T>
T subtract(T, T);

template <class T>
T gen_rot(gsl_rng *r);

template <class R_t, class X_t>
class state_t {
  public:
    D_t D;
    L_t L;
    v_t nv;
    v_t ne;
    graph_t *g;
    double T;
    X_t *spins;
    R_t R;
    double E;
    X_t M; // the "sum" of the spins, like the total magnetization

    std::function <double(X_t, X_t)> J;
    std::function <double(X_t)> H;

    state_t(D_t D, L_t L, double T, std::function <double(X_t, X_t)> J, std::function <double(X_t)> H) : D(D), L(L), T(T), J(J), H(H) {
      graph_t *h = graph_create_square(D, L);
      nv = h->nv;
      ne = h->ne;
      g = graph_add_ext(h);
      graph_free(h);
      spins = (X_t *)malloc(nv * sizeof(X_t));
      for (v_t i = 0; i < nv; i++) {
        init (&(spins[i]));
      }
      init (&R);
      E = - (double)ne * J(spins[0], spins[0]) - (double)nv * H(spins[0]);
      M = scalar_multiple (nv, spins[0]);
    }

    ~state_t() {
      graph_free(g);
      for (v_t i = 0; i < nv; i++) {
        free_spin(spins[i]);
      }
      free(spins);
      free_spin(R);
      free_spin(M);
    }
};

template <q_t q, class T>
struct vector_t { T *x; };

template <q_t q, class T>
void init(vector_t <q, T> *ptr) {
  ptr->x = (T *)calloc(q, sizeof(T));

  ptr->x[0] = (T)1;
}

template <q_t q, class T>
vector_t <q, T> copy (vector_t <q, T> v) {
  vector_t <q, T> v_copy;
 
 v_copy.x = (T *)calloc(q, sizeof(T));

  for (q_t i = 0; i < q; i++) {
    v_copy.x[i] = v.x[i];
  }

  return v_copy;
}

template <q_t q, class T>
void add (vector_t <q, T> v1, vector_t <q, T> v2) {
  for (q_t i = 0; i < q; i++) {
    v1.x[i] += v2.x[i];
  }
}

template <q_t q, class T>
void subtract (vector_t <q, T> v1, vector_t <q, T> v2) {
  for (q_t i = 0; i < q; i++) {
    v1.x[i] -= v2.x[i];
  }
}

template <q_t q, class T>
vector_t <q, T> scalar_multiple(v_t a, vector_t <q, T> v) {
  vector_t <q, T> multiple;
  multiple.x = (T *)malloc(q * sizeof(T));
  for (q_t i = 0; i < q; i++) {
    multiple.x[i] = a * v.x[i];
  }

  return multiple;
}

template <q_t q, class T>
T dot(vector_t <q, T> v1, vector_t <q, T> v2) {
  T prod = 0;

  for (q_t i = 0; i < q; i++) {
    prod += v1.x[i] * v2.x[i];
  }

  return prod;
}

template <q_t q, class T>
void free_spin (vector_t <q, T> v) {
  free(v.x);
}

template <q_t q, class T>
struct orthogonal_t { T *x; };

template <q_t q, class T>
void init(orthogonal_t <q, T> *ptr) {
  ptr->x = (T *)calloc(q * q, sizeof(T));

  for (q_t i = 0; i < q; i++) {
    ptr->x[q * i + i] = (T)1;
  }
}

template <q_t q, class T>
orthogonal_t <q, T> copy (orthogonal_t <q, T> m) {
  orthogonal_t <q, T> m_copy;
  m_copy.x = (T *)calloc(q * q, sizeof(T));

  for (q_t i = 0; i < q * q; i++) {
    m_copy.x[i] = m.x[i];
  }

  return m_copy;
}

template <q_t q, class T>
void free_spin (orthogonal_t <q, T> m) {
  free(m.x);
}

template <q_t q, class T>
vector_t <q, T> act (orthogonal_t <q, T> m, vector_t <q, T> v) {
  vector_t <q, T> v_rot;
  v_rot.x = (T *)calloc(q, sizeof(T));

  for (q_t i = 0; i < q; i++) {
    for (q_t j = 0; j < q; j++) {
      v_rot.x[i] += m.x[q * i + j] * v.x[j];
    }
  }

  return v_rot;
}


template <q_t q, class T>
orthogonal_t <q, T> act (orthogonal_t <q, T> m1, orthogonal_t <q, T> m2) {
  orthogonal_t <q, T> m2_rot;
  m2_rot.x = (T *)calloc(q * q, sizeof(T));

  for (q_t i = 0; i < q; i++) {
    for (q_t j = 0; j < q; j++) {
      for (q_t k = 0; k < q; k++) {
        m2_rot.x[i * q + j] += m1.x[i * q + j] * m2.x[j * q + k];
      }
    }
  }

  return m2_rot;
}

template <q_t q, class T>
vector_t <q, T> act_inverse (orthogonal_t <q, T> m, vector_t <q, T> v) {
  vector_t <q, T> v_rot;
  v_rot.x = (T *)calloc(q, sizeof(T));

  for (q_t i = 0; i < q; i++) {
    for (q_t j = 0; j < q; j++) {
      v_rot.x[i] += m.x[q * j + i] * v.x[j];
    }
  }

  return v_rot;
}

template <q_t q, class T>
orthogonal_t <q, T> act_inverse (orthogonal_t <q, T> m1, orthogonal_t <q, T> m2) {
  orthogonal_t <q, T> m2_rot;
  m2_rot.x = (T *)calloc(q * q, sizeof(T));

  for (q_t i = 0; i < q; i++) {
    for (q_t j = 0; j < q; j++) {
      for (q_t k = 0; k < q; k++) {
        m2_rot.x[i * q + j] += m1.x[j * q + i] * m2.x[j * q + k];
      }
    }
  }

  return m2_rot;
}

template <q_t q>
void generate_rotation (gsl_rng *r, orthogonal_t <q, double> *ptr) {
  double *v = (double *)malloc(q * sizeof(double));
  double v2 = 0;

  for (q_t i = 0; i < q; i++) {
    v[i] = gsl_ran_ugaussian(r);
    v2 += v[i] * v[i];
  }

  ptr->x = (double *)calloc(q * q, sizeof(double));
  
  for (q_t i = 0; i < q; i++) {
    ptr->x[q * i + i] = 1.0;
    for (q_t j = 0; j < q; j++) {
      ptr->x[q * i + j] -= 2 * v[i] * v[j] / v2;
    }
  }

  free(v);
}

template <class R_t, class X_t>
v_t flip_cluster(state_t <R_t, X_t> *state, v_t v0, R_t r, gsl_rng *rand) {
  v_t nv = 0;

  ll_t *stack = NULL;     // create a new stack
  stack_push(&stack, v0); // push the initial vertex to the stack

  bool *marks = (bool *)calloc(state->g->nv, sizeof(bool));

  while (stack != NULL) {
    v_t v = stack_pop(&stack);

    if (!marks[v]) {
      X_t si_old, si_new;
      R_t R_old, R_new;

      si_old = state->spins[v];
      R_old = state->R;

      marks[v] = true;

      if (v == state->g->nv - 1) {
        R_new = act (r, R_old);
      } else {
        si_new = act (r, si_old);
      }

      v_t nn = state->g->v_i[v + 1] - state->g->v_i[v];

      for (v_t i = 0; i < nn; i++) {
        v_t vn = state->g->v_adj[state->g->v_i[v] + i];

        X_t sj;

        if (vn != state->g->nv - 1) {
          sj = state->spins[vn];
        }

        double prob;

        bool is_ext = (v == state->g->nv - 1 || vn == state->g->nv - 1);

        if (is_ext) {
          X_t rs_old, rs_new;
          if (vn == state->g->nv - 1) {
            rs_old = act_inverse (R_old, si_old);
            rs_new = act_inverse (R_old, si_new);
          } else {
            rs_old = act_inverse (R_old, sj);
            rs_new = act_inverse (R_new, sj);
          }
          double dE = state->H(rs_old) - state->H(rs_new);
          prob = 1.0 - exp(-dE / state->T);

          subtract (state->M, rs_old);
          add (state->M, rs_new);
          state->E += dE;

          free_spin (rs_old);
          free_spin (rs_new);
        } else {
          double dE = state->J(si_old, sj) - state->J(si_new, sj);
          prob = 1.0 - exp(-dE / state->T);
          state->E += dE;
        }

        if (gsl_rng_uniform(rand) < prob) { // and with probability...
          stack_push(&stack, vn); // push the neighboring vertex to the stack
        }
      }

      if (v == state->g->nv - 1) {
        free_spin(state->R);
        state->R = R_new;
      } else {
        free_spin(state->spins[v]);
        state->spins[v] = si_new;
      }

      if (v != state->g->nv - 1) { // count the number of non-external sites that flip
        nv++;
      }
    }
  }

  free(marks);

  return nv;
}

