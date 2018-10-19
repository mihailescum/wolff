
#ifndef WOLFF_STATE_H
#define WOLFF_STATE_H

#include <functional>
#include <vector>
#include <random>

#include "graph.hpp"

namespace wolff {

#include "types.h"

template <class R_t, class X_t>
class measurement;

template <class R_t, class X_t>
class system {
  public:
    v_t nv; // number of vertices
    v_t ne; // number of edges
    graph G; // the graph defining the lattice with ghost
    double T; // the temperature
    std::vector<X_t> s; // the state of the ordinary spins

#ifdef WOLFF_BOND_DEPENDENCE
    std::function <double(v_t, const X_t&, v_t, const X_t&)> Z; // coupling between sites
#else
    std::function <double(const X_t&, const X_t&)> Z; // coupling between sites
#endif

#ifndef WOLFF_NO_FIELD
    R_t s0; // the current state of the ghost site
#ifdef WOLFF_SITE_DEPENDENCE
    std::function <double(v_t, const X_t&)> B; // coupling with the external field
#else
    std::function <double(const X_t&)> B; // coupling with the external field
#endif
#endif

    system(graph g, double T,
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
        ) : G(g), T(T), Z(Z)
#ifndef WOLFF_NO_FIELD
             , s0(), B(B)
#endif
    {
      nv = G.nv;
      ne = G.ne;
      s.resize(nv);
#ifndef WOLFF_NO_FIELD
      G.add_ghost();
#endif
#ifdef WOLFF_FINITE_STATES
      finite_states_init(*this);
#endif
    }

    void flip_cluster(v_t, const R_t&, std::mt19937&, measurement<R_t, X_t>&);
    void run_wolff(N_t, std::function <R_t(std::mt19937&, const system<R_t, X_t>&, v_t)> r_gen, measurement<R_t, X_t>& A, std::mt19937& rng);
};

}

#endif

