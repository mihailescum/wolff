
#include "graph.h"

graph_t *graph_create_square(D_t D, L_t L) {
  v_t nv = pow(L, D);
  v_t ne = D * nv;

  v_t *v_i = (v_t *)malloc((nv + 1) * sizeof(v_t));

  for (v_t i = 0; i < nv + 1; i++) {
    v_i[i] = 2 * D * i;
  }

  v_t *v_adj = (v_t *)malloc(2 * D * nv * sizeof(v_t));

  for (v_t i = 0; i < nv; i++) {
    for (D_t j = 0; j < D; j++) {
      v_adj[v_i[i] + 2 * j] = pow(L, j + 1) * (i / ((v_t)pow(L, j + 1))) + fmod(i + pow(L, j), pow(L, j + 1));
      v_adj[v_i[i] + 2 * j + 1] = pow(L, j + 1) * (i / ((v_t)pow(L, j + 1))) + fmod(pow(L, j+1) + i - pow(L, j), pow(L, j + 1));
    }
  }

  graph_t *g = (graph_t *)malloc(sizeof(graph_t));

  g->ne = ne;
  g->nv = nv;
  g->v_i = v_i;
  g->v_adj = v_adj;

  return g;
}

graph_t *graph_add_ext(const graph_t *G) {
  graph_t *tG = (graph_t *)calloc(1, sizeof(graph_t));

  tG->nv = G->nv + 1;
  tG->ne = G->ne + G->nv;

  tG->v_i = (v_t *)malloc((tG->nv + 1) * sizeof(v_t));
  tG->v_adj = (v_t *)malloc(2 * tG->ne * sizeof(v_t));

  for (v_t i = 0; i < G->nv + 1; i++) {
    tG->v_i[i] = G->v_i[i] + i;
  }

  tG->v_i[tG->nv] = 2 * tG->ne;

  for (v_t i = 0; i < G->nv; i++) {
    v_t nn = G->v_i[i + 1] - G->v_i[i];

    for (v_t j = 0; j < nn; j++) {
      tG->v_adj[tG->v_i[i] + j] = G->v_adj[G->v_i[i] + j];
    }

    tG->v_adj[tG->v_i[i] + nn] = G->nv;
    tG->v_adj[tG->v_i[G->nv] + i] = i;
  }

  return tG;
}

void graph_free(graph_t *g) {
  free(g->v_i);
  free(g->v_adj);
  free(g);
}

