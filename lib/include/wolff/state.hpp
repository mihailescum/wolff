
#pragma once

#include <functional>
#include <vector>

#include "types.h"
#include "graph.hpp"

template <class R_t, class X_t>
class state_t {
  private:
    // updating fourier terms F requires many cos and sin calls, faster to do it beforehand.
    std::vector<std::vector<double>> precomputed_cos;
    std::vector<std::vector<double>> precomputed_sin;
  public:
    D_t D;
    L_t L;
    v_t nv; // the number of vertices in the lattice
    v_t ne; // the number of edges in the lattice
    graph_t g; // the graph defining the lattice without ghost
    double T; // the temperature
    std::vector<X_t> spins; // the state of the ordinary spins
#ifndef NOFIELD
    R_t R; // the current state of the ghost site
#endif
    double E; // the system's total energy
    typename X_t::M_t M; // the "sum" of the spins, like the total magnetization
    v_t last_cluster_size; // the size of the last cluster
    std::vector<typename X_t::F_t> ReF;
    std::vector<typename X_t::F_t> ImF;

#ifdef BOND_DEPENDENCE
    std::function <double(v_t, const X_t&, v_t, const X_t&)> J; // coupling between sites
#else
    std::function <double(const X_t&, const X_t&)> J; // coupling between sites
#endif

#ifndef NOFIELD
#ifdef SITE_DEPENDENCE
    std::function <double(v_t, const X_t&)> H; // coupling with the external field
#else
    std::function <double(const X_t&)> H; // coupling with the external field
#endif
#endif

    state_t(D_t D, L_t L, double T,
#ifdef BOND_DEPENDENCE
        std::function <double(v_t, const X_t&, v_t, const X_t&)> J
#else
        std::function <double(const X_t&, const X_t&)> J
#endif
#ifndef NOFIELD
#ifdef SITE_DEPENDENCE
        , std::function <double(v_t, const X_t&)> H
#else
        , std::function <double(const X_t&)> H
#endif
#endif
           ) : D(D), L(L), g(D, L), T(T),
#ifndef NOFIELD
              R(),
#endif
              J(J)
#ifndef NOFIELD
               , H(H)
#endif
    {
      nv = g.nv;
      ne = g.ne;
      spins.resize(nv);
#ifdef BOND_DEPENDENCE
      E = 0;
      for (v_t v = 0; v < nv; v++) {
        for (const v_t &vn : g.v_adj[v]) {
          if (v < vn) {
            E -= J(v, spins[v], vn, spins[vn]);
          }
        }
      }
#else
      E = - (double)ne * J(spins[0], spins[0]);
#endif

#ifndef NOFIELD
      g.add_ext();
#ifdef SITE_DEPENDENCE
      for (v_t i = 0; i < nv; i++) {
        E -= H(i, spins[i]);
      }
#else
      E -= (double)nv * H(spins[0]);
#endif
#endif

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



