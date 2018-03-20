
#include <stdbool.h>
#include <stdlib.h>

#include "types.h"

typedef struct {
  h_t i;
  bool r;
} dihinf_t;

dihinf_t *dihinf_compose(h_t gti, const dihinf_t *g2);

h_t dihinf_act(h_t gi, h_t s);

h_t dihinf_inverse_act(const dihinf_t *g, h_t s);

