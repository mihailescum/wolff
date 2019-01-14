
#ifndef WOLFF_H
#define WOLFF_H

#include <functional>
#include <vector>
#include <random>
#include <cmath>
#include <iterator>
#include <list>
#include <tuple>
#include <queue>

namespace wolff{

  template <class vertex_prop = std::tuple<>, class edge_prop = std::tuple<>>
  class graph {
    /* the graph class describes a lattice on which wolff is run. for most
     * purposes, using the default square lattice constructor is sufficient,
     * but arbitrary lattices can be constructed
     */
    public:
      unsigned D;  // dimension of space
      unsigned L;  // linear size
      unsigned ne; // number of edges
      unsigned nv; // number of vertices

      struct _vertex;

      typedef struct _halfedge {
        struct _vertex &self;     // reference to the vertex this edge comes from
        struct _vertex &neighbor; // reference to the vertex this edge goes to
        edge_prop prop;           // optional properties

        _halfedge(struct _vertex &v1, struct _vertex &v2) : self(v1), neighbor(v2) {};
      } halfedge;

      typedef struct _vertex {
        unsigned ind;              // index of the vertex
        std::list<halfedge> edges; // list of edges incident on this vertex
        vertex_prop prop;          // optional properties
      } vertex;

      std::vector<vertex> vertices; /* vector of vertices, length nv, with 
                                     * vertices[i].ind = i for all 0 <= i < nv
                                     */
      graph() {
        // default constructor for empty graph. use it to build your own!
        D = 0;
        L = 0;
        nv = 0;
        ne = 0;
      };

      graph(unsigned D, unsigned L) : D(D), L(L) {
        // default constructor for square lattice graph
        nv = pow(L, D);
        ne = D * nv;

        vertices.resize(nv);

        for (unsigned i = 0; i < nv; i++) {
          vertices[i].ind = i;
          for (unsigned j = 0; j < D; j++) {
            unsigned n1 = pow(L, j + 1) * (i / ((unsigned)pow(L, j + 1))) +
                                              fmod(i + pow(L, j), pow(L, j + 1));
            unsigned n2 = pow(L, j + 1) * (i / ((unsigned)pow(L, j + 1))) +
                              fmod(pow(L, j + 1) + i - pow(L, j), pow(L, j + 1));

            halfedge f(vertices[i], vertices[n1]);
            halfedge b(vertices[i], vertices[n2]);

            vertices[i].edges.push_back(f);
            vertices[i].edges.push_back(b);
          }
        }
      };

      void add_ghost() {
        // adds a ghost site to any graph
        vertices.resize(nv + 1);

        for (auto it = vertices.begin(); it != std::prev(vertices.end()); it++) {
          halfedge e1(*it, vertices[nv]);
          it->edges.push_back(e1);

          halfedge e2(vertices[nv], *it);
          vertices[nv].edges.push_back(e2);
        }

        vertices[nv].ind = nv;

        ne += nv;
        nv++;
      };
  };

  template <class R_t, class X_t, class G_t = graph<>>
  class measurement;

  template <class R_t, class X_t, class G_t = graph<>>
  class system {
    public:
      unsigned nv; // number of vertices
      unsigned ne; // number of edges
      G_t G; // the graph defining the lattice with ghost
      double T; // the temperature
      std::vector<X_t> s; // the state of the ordinary spins

#ifdef WOLFF_BOND_DEPENDENCE
      std::function <double(const typename G_t::halfedge&, const X_t&, const X_t&)> Z; // coupling between sites
#else
      std::function <double(const X_t&, const X_t&)> Z; // coupling between sites
#endif

#ifndef WOLFF_NO_FIELD
      R_t s0; // the current state of the ghost site
#ifdef WOLFF_SITE_DEPENDENCE
      std::function <double(const typename G_t::vertex&, const X_t&)> B; // coupling with the external field
#else
      std::function <double(const X_t&)> B; // coupling with the external field
#endif
#endif

#ifdef WOLFF_USE_FINITE_STATES
      std::array<std::array<std::array<double, WOLFF_FINITE_STATES_N>, WOLFF_FINITE_STATES_N>, WOLFF_FINITE_STATES_N> Zp;
#ifndef WOLFF_NO_FIELD
      std::array<std::array<double, WOLFF_FINITE_STATES_N>, WOLFF_FINITE_STATES_N> Bp;
#endif
#endif

      system(G_t g, double T,
#ifdef WOLFF_BOND_DEPENDENCE
          std::function <double(const typename G_t::halfedge&, const X_t&, const X_t&)> Z
#else
          std::function <double(const X_t&, const X_t&)> Z
#endif
#ifndef WOLFF_NO_FIELD
#ifdef WOLFF_SITE_DEPENDENCE
          , std::function <double(const typename G_t::vertex&, const X_t&)> B
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
#ifdef WOLFF_USE_FINITE_STATES
        this->finite_states_init();
#endif
      }

      void flip_cluster(typename G_t::vertex& v0, const R_t& r, std::mt19937& rng,
                          measurement<R_t, X_t, G_t>& A) {
        std::uniform_real_distribution<double> dist(0.0, 1.0);

        std::queue<unsigned> queue;
        queue.push(v0.ind);

        std::vector<bool> visited(G.nv, false);

        while (!queue.empty()) {
          unsigned i = queue.front();
          queue.pop();

          if (!visited[i]) { // don't reprocess anyone we've already visited!
            visited[i] = true;

            X_t si_new;
#ifndef WOLFF_NO_FIELD
            R_t s0_new;

            bool we_are_ghost = (i == nv);

            if (we_are_ghost) {
              s0_new = r.act(s0);
            } else
#endif
            {
              si_new = r.act(s[i]);
            }

            for (const typename G_t::halfedge &e : G.vertices[i].edges) {
              double dE, p;
              unsigned j = e.neighbor.ind;

#ifndef WOLFF_NO_FIELD
              bool neighbor_is_ghost = (j == nv);

              if (we_are_ghost || neighbor_is_ghost) {
                X_t s0s_old, s0s_new;
                unsigned non_ghost;

                if (neighbor_is_ghost) {
                  non_ghost = i;
                  s0s_old = s0.act_inverse(s[i]);
                  s0s_new = s0.act_inverse(si_new);
                } else {
                  non_ghost = j;
                  s0s_old = s0.act_inverse(s[j]);
                  s0s_new = s0_new.act_inverse(s[j]);
                }

#ifdef WOLFF_SITE_DEPENDENCE
                dE = B(G.vertices[non_ghost], s0s_old) - B(G.vertices[non_ghost], s0s_new);
#else
                dE = B(s0s_old) - B(s0s_new);
#endif

#ifdef WOLFF_USE_FINITE_STATES
                p = Bp[s0s_old.enumerate()][s0s_new.enumerate()];
#endif

                // run measurement hooks for encountering a ghost bond
                A.ghost_bond_visited(*this, G.vertices[non_ghost], s0s_old, s0s_new, dE);
              } else // this is a perfectly normal bond!
#endif
              {
#ifdef WOLFF_BOND_DEPENDENCE
                dE = Z(e, s[i], s[j]) - Z(e, si_new, s[j]);
#else
                dE = Z(s[i], s[j]) - Z(si_new, s[j]);
#endif

#ifdef WOLFF_USE_FINITE_STATES
                p = Zp[s[i].enumerate()][si_new.enumerate()][s[j].enumerate()];
#endif

                // run measurement hooks for encountering a plain bond
                A.plain_bond_visited(*this, e, si_new, dE);
              }

#ifndef WOLFF_USE_FINITE_STATES
              p = 1.0 - exp(-dE / T);
#endif

              if (dist(rng) < p) {
                queue.push(j); // push the neighboring vertex to the queue 
              }
            }

#ifndef WOLFF_NO_FIELD
            if (we_are_ghost) {
              A.ghost_site_transformed(*this, s0_new);
              s0 = s0_new;
            } else
#endif
            {
              A.plain_site_transformed(*this, G.vertices[i], si_new);
              s[i] = si_new;
            }
          }
        }
      }

      void run_wolff(unsigned N,
          std::function <R_t(std::mt19937&, const system<R_t, X_t, G_t>&, const typename G_t::vertex&)> r_gen, measurement<R_t, X_t, G_t>& A, std::mt19937& rng) {
        std::uniform_int_distribution<unsigned> dist(0, nv - 1);

        for (unsigned n = 0; n < N; n++) {
          unsigned i0 = dist(rng);
          R_t r = r_gen(rng, *this, G.vertices[i0]);

          A.pre_cluster(n, N, *this, G.vertices[i0], r);

          this->flip_cluster(G.vertices[i0], r, rng, A);

          A.post_cluster(n, N, *this);
        }
      }

#ifdef WOLFF_USE_FINITE_STATES
      void finite_states_init() {
#ifndef WOLFF_NO_FIELD
        for (unsigned i = 0; i < WOLFF_FINITE_STATES_N; i++) {
          for (unsigned j = 0; j < WOLFF_FINITE_STATES_N; j++) {
            Bp[i][j] = 1.0 - exp(-(B(X_t(i)) - B(X_t(j))) / T);
          }
        }
#endif
        for (unsigned i = 0; i < WOLFF_FINITE_STATES_N; i++) {
          for (unsigned j = 0; j < WOLFF_FINITE_STATES_N; j++) {
            for (unsigned k = 0; k < WOLFF_FINITE_STATES_N; k++) {
              Zp[i][j][k] = 1.0 - exp(-(Z(X_t(i), X_t(k)) - Z(X_t(j), X_t(k))) / T);
            }
          }
        }
      }
#endif

  };

  template <class R_t, class X_t, class G_t>
  class measurement {
    public:
      virtual void pre_cluster(unsigned, unsigned, const system<R_t, X_t, G_t>&, const typename G_t::vertex& v, const R_t&) {};

      virtual void plain_bond_visited(const system<R_t, X_t, G_t>&, const typename G_t::halfedge& e, const X_t&, double) {};
      virtual void plain_site_transformed(const system<R_t, X_t, G_t>&, const typename G_t::vertex& v, const X_t&) {};

#ifndef WOLFF_NO_FIELD
      virtual void ghost_bond_visited(const system<R_t, X_t, G_t>&, const typename G_t::vertex& v, const X_t&, const X_t&, double) {};
      virtual void ghost_site_transformed(const system<R_t, X_t, G_t>&, const R_t&) {};
#endif

      virtual void post_cluster(unsigned, unsigned, const system<R_t, X_t, G_t>&) {};
  };

}

#endif

