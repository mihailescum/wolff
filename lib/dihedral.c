
#include "dihedral.h"

dihedral_t *dihedral_compose(q_t q, q_t g1i, const dihedral_t *g2) {
  // we only need to consider the action of reflections
  dihedral_t *g3 = (dihedral_t *)malloc(1 * sizeof(dihedral_t));

  g3->r = !g2->r;
  g3->i = (g1i + q - g2->i) % q;

  return g3;
}

q_t dihedral_act(q_t q, q_t gi, q_t s) {
  // we only need to consider the action of reflections

  return (gi + q - s) % q;
}

q_t dihedral_inverse_act(q_t q, const dihedral_t *g, q_t s) {
  if (g->r) {
    return (q - ((q + s - g->i) % q)) % q;
  } else {
    return (q + s - g->i) % q;
  }
}


