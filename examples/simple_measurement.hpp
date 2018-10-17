
#include <wolff/measurement.hpp>

template <class R_t, class X_t>
class simple_measurement : public wolff_measurement<R_t, X_t> {
  private:
    N_t n;

    double E;
    typename X_t::M_t M;
    v_t C;

    double totalE;
    typename X_t::F_t totalM;
    double totalC;

  public:
    simple_measurement(const wolff_system<R_t, X_t>& S) {
      n = 0;
      M = S.nv * S.s[0];
      E = - (S.ne * S.Z(S.s[0], S.s[0]) + S.nv * S.B(S.s[0]));

      totalE = 0;
      totalM = 0.0 * (S.s[0]);
      totalC = 0;
    }

    void pre_cluster(N_t, N_t, const wolff_system<R_t, X_t>&, v_t, const R_t&) {
      C = 0;
    }

    void plain_bond_visited(const wolff_system<R_t, X_t>&, v_t, const X_t&, v_t, double dE) {
      E += dE;
    }

    void ghost_bond_visited(const wolff_system<R_t, X_t>&, v_t, const X_t& s_old, const X_t& s_new, double dE) {
      E += dE;
      M += s_new - s_old;
    }

    void plain_site_transformed(const wolff_system<R_t, X_t>&, v_t,  const X_t&) {
      C++;
    }

    void ghost_site_transformed(const wolff_system<R_t, X_t>&, const R_t&) {}

    void post_cluster(N_t, N_t, const wolff_system<R_t, X_t>&) {
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

