
#pragma once

#include <wolff/state.hpp>
#include <wolff/meas.h>
#include "correlation.hpp"
#include <functional>

#define POSSIBLE_MEASUREMENTS 4
const unsigned char measurement_energy        = 1 << 0;
const unsigned char measurement_clusterSize   = 1 << 1;
const unsigned char measurement_magnetization = 1 << 2;
const unsigned char measurement_fourierZero    = 1 << 3;

char const *measurement_labels[] = {"E", "S", "M", "F"};

template <class R_t, class X_t>
class wolff_research_measurements : public wolff_measurement<R_t, X_t> {
  private:
    std::vector<std::vector<double>> precomputed_sin;
    std::vector<std::vector<double>> precomputed_cos;
    FILE **files;
    D_t D;
    unsigned char flags;
    std::function <void(const state_t <R_t, X_t>&, const wolff_research_measurements<R_t, X_t>&)> other_f;
    bool silent;
  public:
    v_t last_cluster_size;
    double E;
    typename X_t::M_t M;
    std::vector<typename X_t::F_t> ReF;
    std::vector<typename X_t::F_t> ImF;

    wolff_research_measurements(unsigned char flags, unsigned long timestamp, std::function <void(const state_t <R_t, X_t>&, const wolff_research_measurements<R_t, X_t>&)> other_f, state_t<R_t, X_t>& s, bool silent) : flags(flags), D(s.D), other_f(other_f), silent(silent) {
      FILE **files = (FILE **)calloc(POSSIBLE_MEASUREMENTS, sizeof(FILE *));

      for (uint8_t i = 0; i < POSSIBLE_MEASUREMENTS; i++) {
        if (flags & (1 << i)) {
          char *filename = (char *)malloc(255 * sizeof(char));
          sprintf(filename, "wolff_%lu_%s.dat", timestamp, measurement_labels[i]);
          files[i] = fopen(filename, "wb");
          free(filename);
        }
      }

#ifdef BOND_DEPENDENCE
      E = 0;
      for (v_t v = 0; v < s.nv; v++) {
        for (const v_t &vn : s.g.v_adj[v]) {
          if (v < vn) {
            E -= s.J(v, s.spins[v], vn, s.spins[vn]);
          }
        }
      }
#else
      E = - (double)s.ne * s.J(s.spins[0], s.spins[0]);
#endif

#ifndef NOFIELD
#ifdef SITE_DEPENDENCE
      for (v_t i = 0; i < s.nv; i++) {
        E -= s.H(i, s.spins[i]);
      }
#else
      E -= (double)s.nv * s.H(s.spins[0]);
#endif
#endif

      M = s.spins[0] * s.nv;

      ReF.resize(s.D);
      ImF.resize(s.D);
      for (D_t i = 0; i < s.D; i++) {
        ReF[i] = s.spins[0] * 0.0;
        ImF[i] = s.spins[0] * 0.0;
      }
      precomputed_cos.resize(s.nv);
      precomputed_sin.resize(s.nv);
      for (v_t i = 0; i < s.nv; i++) {
        precomputed_cos[i].resize(s.D);
        precomputed_sin[i].resize(s.D);
        for (D_t j = 0; j < s.D; j++) {
          precomputed_cos[i][j] = cos(2 * M_PI * s.g.coordinate[i][j] / (double)s.L);
          precomputed_sin[i][j] = sin(2 * M_PI * s.g.coordinate[i][j] / (double)s.L);
        }
      }

      if (!silent) printf("\n");

    }

    void pre_cluster(const state_t<R_t, X_t>& s, count_t step, count_t N, v_t v0, const R_t& R) {
      last_cluster_size = 0;
      if (!silent) printf("\033[F\033[JWOLFF: step %" PRIu64 " / %" PRIu64 ": E = %.2f, S = %" PRIv "\n", step, N, E, last_cluster_size);

    }

    void plain_bond_added(v_t vi, const X_t& si_old, const X_t& si_new, v_t vn, const X_t& sn, double dE) {
      E += dE;
    }

    void ghost_bond_added(v_t v, const X_t& rs_old, const X_t& rs_new, double dE) {
      E += dE;
      M += rs_new - rs_old;

#ifdef DIMENSION
      for (D_t i = 0; i < DIMENSION; i++)
#else
      for (D_t i = 0; i < D; i++)
#endif
      {
        ReF[i] += (rs_new - rs_old) * precomputed_cos[v][i];
        ImF[i] += (rs_new - rs_old) * precomputed_sin[v][i];
      }
    }

    void plain_site_transformed(v_t v, const X_t& s_old, const X_t& s_new) {
      last_cluster_size++;
    }

    void ghost_site_transformed(const R_t& r_old, const R_t& r_new) {
    }

    void post_cluster(const state_t<R_t, X_t>& s, count_t step, count_t N) {
      if (flags & measurement_energy) {
        float smaller_E = (float)E;
        fwrite(&smaller_E, sizeof(float), 1, files[0]);
      }
      if (flags & measurement_clusterSize) {
        fwrite(&(last_cluster_size), sizeof(uint32_t), 1, files[1]);
      }
      if (flags & measurement_magnetization) {
        write_magnetization(M, files[2]);
      }
      if (flags & measurement_fourierZero) {
        float smaller_X = (float)correlation_length<X_t>(ReF, ImF, D);
        fwrite(&smaller_X, sizeof(float), 1, files[3]);
      }

      other_f(s, *this);

    }

    ~wolff_research_measurements() {
      for (uint8_t i = 0; i < POSSIBLE_MEASUREMENTS; i++) {
        if (flags & (1 << i)) {
          fclose(files[i]);
        }
      }

      if (!silent) {
        printf("\033[F\033[J");
      }
      printf("WOLFF COMPLETE: E = %.2f, S = %" PRIv "\n", E, last_cluster_size);

      free(files);
    }
};

