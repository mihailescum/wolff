
#include "dihedral.h"

dihedral_t *dihedral_compose(q_t q, q_t g1i, const dihedral_t *g2) {
  // we only need to consider the action of reflections
  dihedral_t *g3 = (dihedral_t *)malloc(1 * sizeof(dihedral_t));

  g3->r = !g2->r;
  g3->i = (g1i + q - g2->i) % q;

  return g3;
}

q_t dihedral_act(q_t q, q_t gi, bool r, q_t s) {
  // we only need to consider the action of reflections

  if (r) {
    return (gi + q - s) % q;
  } else {
    return (gi + s) % q;
  }
}

q_t dihedral_inverse_act(q_t q, const dihedral_t *g, q_t s) {
  if (g->r) {
    return (q - ((q + s - g->i) % q)) % q;
  } else {
    return (q + s - g->i) % q;
  }
}

q_t *dihedral_gen_transformations(q_t q) {
  q_t *transformations = (q_t *)malloc(2 * q * q * sizeof(q_t));

  for (q_t i = 0; i < q; i++) {
    for (q_t j = 0; j < q; j++) {
      transformations[q * i + j] = dihedral_act(q, i, false, j);
      transformations[q * q + q * i + j] = dihedral_act(q, i, true, j);
    }
  }

  return transformations;
}

R_t *dihedral_gen_involutions(q_t q) {
  R_t *transformations = (R_t *)malloc(q * sizeof(R_t));

  for (q_t i = 0; i < q; i++) {
    transformations[i] = q + i;
  }

  return transformations;
}


