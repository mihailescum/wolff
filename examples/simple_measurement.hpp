
#include <wolff/measurement.hpp>

using namespace wolff;

template <class R_t, class X_t>
class simple_measurement : public measurement<R_t, X_t> {
  private:
    N_t n;

    double E;
    typename X_t::M_t M;
    v_t C;

    double totalE;
    typename X_t::F_t totalM;
    double totalC;

  public:
    simple_measurement(const system<R_t, X_t>& S) {
      n = 0;
      M = S.nv * S.s[0];
      E = 0;

#ifdef WOLFF_BOND_DEPENDENCE
      for (v_t i = 0; i < S.nv; i++) {
        for (v_t j : S.G.adjacency[i]) {
          E -= 0.5 * S.Z(i, S.s[i], j, S.s[j]);
        }
      }
#else
      E -= S.ne * S.Z(S.s[0], S.s[0]);
#endif

#ifndef WOLFF_NO_FIELD
#ifdef WOLFF_SITE_DEPENDENCE
      for (v_t i = 0; i < S.nv; i++) {
        E -= S.B(i, S.s[i]);
      }
#else
      E -= S.nv * S.B(S.s[0]);
#endif
#endif

      totalE = 0;
      totalM = 0.0 * (S.s[0]);
      totalC = 0;
    }

    void pre_cluster(N_t, N_t, const system<R_t, X_t>&, v_t, const R_t&) {
      C = 0;
    }

    void plain_bond_visited(const system<R_t, X_t>&, v_t, const X_t&, v_t, double dE) {
      E += dE;
    }

    void ghost_bond_visited(const system<R_t, X_t>&, v_t, const X_t& s_old, const X_t& s_new, double dE) {
      E += dE;
      M += s_new - s_old;
    }

    void plain_site_transformed(const system<R_t, X_t>& S, v_t i, const X_t& si_new) {
      C++;

#ifdef WOLFF_NO_FIELD
      M += si_new - S.s[i];
#endif
    }

    void post_cluster(N_t, N_t, const system<R_t, X_t>&) {
      totalE += E;
      totalM += M;
      totalC += C;
      n++;
    }

    double avgE() {
      return totalE / n;
    }

    typename X_t::F_t avgM() {
      return totalM / n;
    }

    double avgC() {
      return totalC / n;
    }
};

