
#include "dihinf.h"

dihinf_t *dihinf_compose(h_t g1i, const dihinf_t *g2) {
  // we only need to consider the action of reflections
  dihinf_t *g3 = (dihinf_t *)malloc(1 * sizeof(dihinf_t));

  g3->r = !g2->r;
  g3->i = g1i - g2->i;

  return g3;
}

h_t dihinf_act(h_t gi, h_t s) {
  // we only need to consider the action of reflections

  return gi - s;
}

h_t dihinf_inverse_act(const dihinf_t *g, h_t s) {
  if (g->r) {
    return g->i - s;
  } else {
    return s - g->i;
  }
}


