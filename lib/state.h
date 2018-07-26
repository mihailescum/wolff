
#pragma once

#include <functional>
#include <vector>

#include "types.h"
#include "graph.h"

template <class R_t, class X_t>
class state_t {
  public:
    D_t D;
    L_t L;
    v_t nv;
    v_t ne;
    graph_t g;
    double T;
    std::vector<X_t> spins;
    R_t R;
    double E;
    typename X_t::M_t M; // the "sum" of the spins, like the total magnetization
    v_t last_cluster_size;
    typename X_t::F_t *ReF;
    typename X_t::F_t *ImF;
    // updating fourier terms F requires many cos and sin calls, faster to do it beforehand.
    double *precomputed_cos;
    double *precomputed_sin;

    std::function <double(X_t, X_t)> J;
    std::function <double(X_t)> H;

    state_t(D_t D, L_t L, double T, std::function <double(X_t, X_t)> J, std::function <double(X_t)> H) : D(D), L(L), g(D, L), T(T), J(J), H(H) {
      nv = g.nv;
      ne = g.ne;
      g.add_ext();
      spins.resize(nv);
      init (&R);
      E = - (double)ne * J(spins[0], spins[0]) - (double)nv * H(spins[0]);
      M = spins[0] * nv;
      last_cluster_size = 0;
      ReF = (typename X_t::F_t *)malloc(D * sizeof(typename X_t::F_t));
      ImF = (typename X_t::F_t *)malloc(D * sizeof(typename X_t::F_t));
      for (D_t i = 0; i < D; i++) {
        ReF[i] = spins[0] * 0.0;
        ImF[i] = spins[0] * 0.0;
      }
      precomputed_cos = (double *)malloc(L * sizeof(double));
      precomputed_sin = (double *)malloc(L * sizeof(double));
      for (L_t i = 0; i < L; i++) {
        precomputed_cos[i] = cos(2 * M_PI * (double)i / (double)L);
        precomputed_sin[i] = sin(2 * M_PI * (double)i / (double)L);
      }
    }

    ~state_t() {
      free_spin(R);
      free(ReF);
      free(ImF);
      free(precomputed_sin);
      free(precomputed_cos);
    }
};



