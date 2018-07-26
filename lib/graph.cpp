
#include "graph.h"

graph_t::graph_t(D_t D, L_t L) {
  nv = pow(L, D);
  ne = D * nv;

  v_adj.resize(nv);

  for (std::vector<v_t> v_adj_i : v_adj) {
    v_adj_i.reserve(2 * D);
  }

  for (v_t i = 0; i < nv; i++) {
    for (D_t j = 0; j < D; j++) {
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
  v_adj[nv].reserve(nv);

  for (v_t i = 0; i < nv; i++) {
    v_adj[nv].push_back(i);
  }

  ne += nv;
  nv += 1;
}

