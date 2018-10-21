
#ifndef WOLFF_MEASUREMENTS_H
#define WOLFF_MEASUREMENTS_H

#include "system.hpp"

namespace wolff {

template <class R_t, class X_t>
class measurement {
  public:
    virtual void pre_cluster(N_t, N_t, const system<R_t, X_t>&, v_t, const R_t&) {};

    virtual void plain_bond_visited(const system<R_t, X_t>&, v_t, const X_t&, v_t, double) {};
    virtual void plain_site_transformed(const system<R_t, X_t>&, v_t, const X_t&) {};

#ifndef WOLFF_NO_FIELD
    virtual void ghost_bond_visited(const system<R_t, X_t>&, v_t, const X_t&, const X_t&, double) {};
    virtual void ghost_site_transformed(const system<R_t, X_t>&, const R_t&) {};
#endif

    virtual void post_cluster(N_t, N_t, const system<R_t, X_t>&) {};
};

}

#endif

