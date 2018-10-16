
#pragma once

#include <functional>
#include <vector>

#include "types.h"
#include "graph.hpp"

template <class R_t, class X_t>
class state_t {
  public:
    D_t D; // the dimension of the system
    L_t L; // the linear size of the lattice
    v_t nv; // the number of vertices in the original lattice
    v_t ne; // the number of edges in the original lattice
    graph_t g; // the graph defining the lattice with ghost
    double T; // the temperature
    std::vector<X_t> spins; // the state of the ordinary spins
#ifndef NOFIELD
    R_t R; // the current state of the ghost site
#endif

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
        , lattice_t lat = SQUARE_LATTICE) : D(D), L(L), g(D, L, lat), T(T),
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
#ifndef NOFIELD
      g.add_ext();
#endif
    }
};

