
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
    std::vector<typename X_t::F_t> ReF;
    std::vector<typename X_t::F_t> ImF;
    // updating fourier terms F requires many cos and sin calls, faster to do it beforehand.
    std::vector<double> precomputed_cos;
    std::vector<double> precomputed_sin;

    std::function <double(const X_t&, const X_t&)> J;
    std::function <double(const X_t&)> H;

    state_t(D_t D, L_t L, double T, std::function <double(const X_t&, const X_t&)> J, std::function <double(const X_t&)> H) : D(D), L(L), g(D, L), T(T), R(), J(J), H(H) {
      nv = g.nv;
      ne = g.ne;
      g.add_ext();
      spins.resize(nv);
      E = - (double)ne * J(spins[0], spins[0]) - (double)nv * H(spins[0]);
      M = spins[0] * nv;
      last_cluster_size = 0;
      ReF.resize(D);
      ImF.resize(D);
      for (D_t i = 0; i < D; i++) {
        ReF[i] = spins[0] * 0.0;
        ImF[i] = spins[0] * 0.0;
      }
      precomputed_cos.resize(L);
      precomputed_sin.resize(L);
      for (L_t i = 0; i < L; i++) {
        precomputed_cos[i] = cos(2 * M_PI * (double)i / (double)L);
        precomputed_sin[i] = sin(2 * M_PI * (double)i / (double)L);
      }
    }
};



