
#include "graph.h"

graph_t::graph_t(D_t D, L_t L) {
  nv = pow(L, D);
  ne = D * nv;

  v_adj.resize(nv);
  coordinate.resize(nv);

  for (std::vector<v_t> v_adj_i : v_adj) {
    v_adj_i.reserve(2 * D);
  }

  for (v_t i = 0; i < nv; i++) {
    coordinate[i].resize(D);
    for (D_t j = 0; j < D; j++) {
      coordinate[i][j] = (i / (v_t)pow(L, D - j - 1)) % L;

      v_adj[i].push_back(pow(L, j + 1) * (i / ((v_t)pow(L, j + 1))) + fmod(i + pow(L, j), pow(L, j + 1)));
      v_adj[i].push_back(pow(L, j + 1) * (i / ((v_t)pow(L, j + 1))) + fmod(pow(L, j+1) + i - pow(L, j), pow(L, j + 1)));
    }
  }
}

void graph_t::add_ext() {
  for (std::vector<v_t>& v_adj_i : v_adj) {
    v_adj_i.push_back(nv);
  }

  v_adj.resize(nv + 1);
  coordinate.resize(nv + 1);
  v_adj[nv].reserve(nv);

  for (v_t i = 0; i < nv; i++) {
    v_adj[nv].push_back(i);
  }

  coordinate[nv].resize(coordinate[0].size());

  ne += nv;
  nv += 1;
}

