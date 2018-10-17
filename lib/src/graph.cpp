
#include <wolff/graph.hpp>

graph_t::graph_t(D_t D, L_t L, lattice_t lat) {
  switch (lat) {
    case SQUARE_LATTICE: {
        nv = pow(L, D);
        ne = D * nv;

        adj.resize(nv);
        coordinate.resize(nv);

        for (std::vector<v_t> adj_i : adj) {
          adj_i.reserve(2 * D);
        }

        for (v_t i = 0; i < nv; i++) {
          coordinate[i].resize(D);
          for (D_t j = 0; j < D; j++) {
            coordinate[i][j] = (i / (v_t)pow(L, D - j - 1)) % L;

            adj[i].push_back(pow(L, j + 1) * (i / ((v_t)pow(L, j + 1))) + fmod(i + pow(L, j), pow(L, j + 1)));
            adj[i].push_back(pow(L, j + 1) * (i / ((v_t)pow(L, j + 1))) + fmod(pow(L, j+1) + i - pow(L, j), pow(L, j + 1)));
          }
        }
        break;
      }
    case DIAGONAL_LATTICE: {
        nv = D * pow(L, D);
        ne = 2 * nv;

        adj.resize(nv);
        coordinate.resize(nv);

        for (std::vector<v_t> adj_i : adj) {
          adj_i.reserve(4 * (D - 1));
        }

        for (D_t i = 0; i < D; i++) {
          v_t sb = i * pow(L, D);

          for (v_t j = 0; j < pow(L, D); j++) {
            v_t vc = sb + j;

            adj[vc].push_back(((i + 1) % D) * pow(L, D) + j);
            adj[vc].push_back(((i + 1) % D) * pow(L, D) + pow(L, D - 1) * (j / (v_t)pow(L, D - 1)) + (j + 1 - 2 * (i % 2)) % L);
            adj[vc].push_back(((i + 1) % D) * pow(L, D) + pow(L, D - 1) * ((L + (j/ (v_t)pow(L, D - 1)) - 1 + 2 * (i % 2)) % L) + (j - i) % L);
            adj[vc].push_back(((i + 1) % D) * pow(L, D) + pow(L, D - 1) * ((L + (j/ (v_t)pow(L, D - 1)) - 1 + 2 * (i % 2)) % L) + (j + 1 - i) % L);
          }
        }
        break;
      }
  }
}

void graph_t::add_ext() {
  for (std::vector<v_t>& adj_i : adj) {
    adj_i.push_back(nv);
  }

  adj.resize(nv + 1);
  coordinate.resize(nv + 1);
  adj[nv].reserve(nv);

  for (v_t i = 0; i < nv; i++) {
    adj[nv].push_back(i);
  }

  coordinate[nv].resize(coordinate[0].size());

  ne += nv;
  nv += 1;
}

