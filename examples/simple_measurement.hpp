
#include <wolff.hpp>

using namespace wolff;

template <class R_t, class X_t, class G_t>
class simple_measurement : public measurement<R_t, X_t, G_t> {
  private:
    unsigned n;

    double E;
    typename X_t::M_t M;
    unsigned C;

    double totalE;
    typename X_t::F_t totalM;
    double totalC;

  public:
    simple_measurement(const wolff::system<R_t, X_t, G_t>& S) {
      n = 0;
      M = S.nv * S.s[0];
      E = 0;

#ifdef WOLFF_BOND_DEPENDENCE
      for (unsigned i = 0; i < S.nv; i++) {
        for (const typename G_t::halfedge& e : S.G.vertices[i].edges) {
          E -= 0.5 * S.Z(e, S.s[i], S.s[j]);
        }
      }
#else
      E -= S.ne * S.Z(S.s[0], S.s[0]);
#endif

#ifndef WOLFF_NO_FIELD
#ifdef WOLFF_SITE_DEPENDENCE
      for (unsigned i = 0; i < S.nv; i++) {
        E -= S.B(S.G.vertices[i], S.s[i]);
      }
#else
      E -= S.nv * S.B(S.s[0]);
#endif
#endif

      totalE = 0;
      totalM = 0.0 * (S.s[0]);
      totalC = 0;
    }

    void pre_cluster(unsigned, unsigned, const wolff::system<R_t, X_t, G_t>&, const typename G_t::vertex&, const R_t&) override {
      C = 0;
    }

    void plain_bond_visited(const wolff::system<R_t, X_t, G_t>&, const typename G_t::halfedge&, const X_t&, double dE) override {
      E += dE;
    }

#ifndef WOLFF_NO_FIELD
    void ghost_bond_visited(const wolff::system<R_t, X_t, G_t>&, const typename G_t::vertex&, const X_t& s_old, const X_t& s_new, double dE) override {
      E += dE;
      M += s_new - s_old;
    }
#endif

    void plain_site_transformed(const wolff::system<R_t, X_t, G_t>& S, const typename G_t::vertex& v, const X_t& si_new) override {
      C++;

#ifdef WOLFF_NO_FIELD
      M += si_new - S.s[v.ind];
#endif
    }

    void post_cluster(unsigned, unsigned, const wolff::system<R_t, X_t, G_t>&) override {
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

