
#include <stdbool.h>
#include <stdlib.h>

#include "types.h"

typedef struct {
  q_t i;
  bool r;
} dihedral_t;

dihedral_t *dihedral_compose(q_t q, q_t gti, const dihedral_t *g2);

q_t dihedral_act(q_t q, q_t gi, q_t s);

q_t dihedral_inverse_act(q_t q, const dihedral_t *g, q_t s);

