
#pragma once

#include "types.h"

template <class R_t, class X_t>
class wolff_measurement {
  public:
    virtual void pre_cluster(const state_t<R_t, X_t>&, count_t, count_t, v_t, const R_t&) = 0;

    virtual void plain_bond_added(v_t, const X_t&, const X_t&, v_t, const X_t&, double) = 0;
    virtual void ghost_bond_added(v_t, const X_t&, const X_t&, double) = 0;

    virtual void plain_site_transformed(v_t, const X_t&, const X_t&) = 0;
    virtual void ghost_site_transformed(const R_t&, const R_t&) = 0;

    virtual void post_cluster(const state_t<R_t, X_t>&, count_t, count_t) = 0;
};

