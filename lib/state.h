
#pragma once

#include <functional>

#include "types.h"
#include "graph.h"

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



