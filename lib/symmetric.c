
#include "symmetric.h"

q_t *symmetric_compose(q_t q, const q_t *g1, const q_t *g2) {
  q_t *g3 = (q_t *)malloc(q * sizeof(q_t));

  for (q_t i = 0; i < q; i++) {
    g3[i] = g1[g2[i]];
  }

  return g3;
}

q_t symmetric_act(const q_t *g, q_t s) {
  return g[s];
}

q_t *symmetric_invert(q_t q, const q_t *g) {
  q_t *g_inv = (q_t *)malloc(q * sizeof(q_t));

  for (q_t i = 0; i < q; i++) {
    g_inv[g[i]] = i;
  }

  return g_inv;
}

