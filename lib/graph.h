
#pragma once

#include <inttypes.h>
#include <math.h>
#include <stdlib.h>

#include "types.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
  v_t ne;
  v_t nv;
  v_t *v_i;
  v_t *v_adj;
} graph_t;

graph_t *graph_create_square(D_t D, L_t L);
graph_t *graph_add_ext(const graph_t *G);
void graph_free(graph_t *h);

#ifdef __cplusplus
}
#endif

