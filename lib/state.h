
#pragma once

#include <functional>
#include <vector>

#include "types.h"
#include "graph.h"

template <class R_t, class X_t>
class state_t {
  private:
    // updating fourier terms F requires many cos and sin calls, faster to do it beforehand.
    std::vector<std::vector<double>> precomputed_cos;
    std::vector<std::vector<double>> precomputed_sin;
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
      precomputed_cos.resize(nv);
      precomputed_sin.resize(nv);
      for (v_t i = 0; i < nv; i++) {
        precomputed_cos[i].resize(D);
        precomputed_sin[i].resize(D);
        for (D_t j = 0; j < D; j++) {
          precomputed_cos[i][j] = cos(2 * M_PI * g.coordinate[i][j] / (double)L);
          precomputed_sin[i][j] = sin(2 * M_PI * g.coordinate[i][j] / (double)L);
        }
      }
    }

    void update_magnetization(const X_t& s_old, const X_t& s_new) {
      M += s_new - s_old;
    }

    void update_energy(const double& dE) {
      E += dE;
    }

    void update_fourierZero(v_t v, const X_t& s_old, const X_t& s_new) {
#ifdef DIMENSION
      for (D_t i = 0; i < DIMENSION; i++) {
#else
      for (D_t i = 0; i < D; i++) {
#endif
        ReF[i] += (s_new - s_old) * precomputed_cos[v][i];
        ImF[i] += (s_new - s_old) * precomputed_sin[v][i];
      }
    }
};



