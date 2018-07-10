
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

void swap(q_t *q1, q_t *q2) {
  q_t temp = *q1;
  *q1 = *q2;
  *q2 = temp;
}

R_t factorial(q_t q) {
  if (q == 0) {
    return 1;
  } else {
    return q * factorial(q - 1);
  }
}

void permute(q_t *a, q_t l, q_t r, R_t pos, q_t *transformations) {
  if (l == r - 1) {
    for (q_t i = 0; i < r; i++) {
      transformations[r * pos + i] = a[i];
    }
  } else {
    for (q_t i = l; i < r; i++) {
      swap((a+l), (a+i));
      permute(a, l+1, r, pos + (i - l) * factorial(r - l - 1), transformations);
      swap((a+l), (a+i));
    }
  }
}

q_t *symmetric_gen_transformations(q_t q) {
  q_t *transformations = (q_t *)malloc(q * factorial(q) * sizeof(q_t));
  q_t *tmp = (q_t *)malloc(q * sizeof(q_t));

  for (q_t i = 0; i < q; i++) {
    tmp[i] = i;
  }

  permute(tmp, 0, q, 0, transformations);

  free(tmp);

  return transformations;
}

