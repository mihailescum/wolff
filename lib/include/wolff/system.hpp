
#ifndef WOLFF_STATE_H
#define WOLFF_STATE_H

#include <functional>
#include <vector>

#include "types.h"
#include "graph.hpp"

template <class R_t, class X_t>
class wolff_system {
  public:
    D_t D; // dimension
    L_t L; // linear size
    v_t nv; // number of vertices
    v_t ne; // number of edges
    graph_t G; // the graph defining the lattice with ghost
    double T; // the temperature
    std::vector<X_t> s; // the state of the ordinary spins
#ifndef WOLFF_NO_FIELD
    R_t s0; // the current state of the ghost site
#endif

#ifdef WOLFF_BOND_DEPENDENCE
    std::function <double(v_t, const X_t&, v_t, const X_t&)> Z; // coupling between sites
#else
    std::function <double(const X_t&, const X_t&)> Z; // coupling between sites
#endif

#ifndef WOLFF_NO_FIELD
#ifdef WOLFF_SITE_DEPENDENCE
    std::function <double(v_t, const X_t&)> B; // coupling with the external field
#else
    std::function <double(const X_t&)> B; // coupling with the external field
#endif
#endif

    wolff_system(D_t D, L_t L, double T,
#ifdef WOLFF_BOND_DEPENDENCE
        std::function <double(v_t, const X_t&, v_t, const X_t&)> Z
#else
        std::function <double(const X_t&, const X_t&)> Z
#endif
#ifndef WOLFF_NO_FIELD
#ifdef WOLFF_SITE_DEPENDENCE
        , std::function <double(v_t, const X_t&)> B
#else
        , std::function <double(const X_t&)> B
#endif
#endif
        , lattice_t lat = SQUARE_LATTICE) : D(D), L(L), G(D, L, lat), T(T), Z(Z)
#ifndef WOLFF_NO_FIELD
             , s0(), B(B)
#endif
    {
      nv = G.nv;
      ne = G.ne;
      s.resize(nv);
#ifndef WOLFF_NO_FIELD
      G.add_ext();
#endif
#ifdef WOLFF_FINITE_STATES
      finite_states_init(*this);
#endif
    }
};

#endif

