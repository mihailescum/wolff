
#ifndef WOLFF_MEASUREMENTS
#define WOLFF_MEASUREMENTS

#include "system.hpp"

template <class R_t, class X_t>
class wolff_measurement {
  public:
    virtual void pre_cluster(N_t, N_t, const wolff_system<R_t, X_t>&, v_t, const R_t&) = 0;

    virtual void plain_bond_visited(const wolff_system<R_t, X_t>&, v_t, const X_t&, v_t, double) = 0;
    virtual void plain_site_transformed(const wolff_system<R_t, X_t>&, v_t, const X_t&) = 0;

#ifndef WOLFF_NO_FIELD
    virtual void ghost_bond_visited(const wolff_system<R_t, X_t>&, v_t, const X_t&, const X_t&, double) = 0;
    virtual void ghost_site_transformed(const wolff_system<R_t, X_t>&, const R_t&) = 0;
#endif

    virtual void post_cluster(N_t, N_t, const wolff_system<R_t, X_t>&) = 0;
};

#endif

