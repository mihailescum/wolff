
#pragma once

#include <stdlib.h>

#include "types.h"

q_t *symmetric_compose(q_t q, const q_t *g1, const q_t *g2);

q_t symmetric_act(const q_t *g, q_t s);

q_t *symmetric_invert(q_t q, const q_t *g);

q_t *symmetric_gen_transformations(q_t q);

