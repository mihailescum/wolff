
#include <wolff/graph.hpp>

namespace wolff {

graph::graph() {
  D = 0;
  L = 0;
  nv = 0;
  ne = 0;
}

graph::graph(D_t D, L_t L) : D(D), L(L) {
  nv = pow(L, D);
  ne = D * nv;

  adjacency.resize(nv);
  coordinate.resize(nv);

  for (std::vector<v_t> adj_i : adjacency) {
    adj_i.reserve(2 * D);
  }

  for (v_t i = 0; i < nv; i++) {
    coordinate[i].resize(D);
    for (D_t j = 0; j < D; j++) {
      coordinate[i][j] = (i / (v_t)pow(L, D - j - 1)) % L;

      adjacency[i].push_back(pow(L, j + 1) * (i / ((v_t)pow(L, j + 1))) + fmod(i + pow(L, j), pow(L, j + 1)));
      adjacency[i].push_back(pow(L, j + 1) * (i / ((v_t)pow(L, j + 1))) + fmod(pow(L, j+1) + i - pow(L, j), pow(L, j + 1)));
    }
  }
}

void graph::add_ghost() {
  for (std::vector<v_t>& adj_i : adjacency) {
    adj_i.push_back(nv);
  }

  adjacency.resize(nv + 1);
  coordinate.resize(nv + 1);
  adjacency[nv].reserve(nv);

  for (v_t i = 0; i < nv; i++) {
    adjacency[nv].push_back(i);
  }

  coordinate[nv].resize(coordinate[0].size());

  ne += nv;
  nv++;
}

}

