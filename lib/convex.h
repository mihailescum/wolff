
#pragma once

#include <inttypes.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

#include "types.h"

typedef struct {
  count_t x;
  double y;
} point_t;

typedef struct list_tag {
  struct list_tag *prev;
  struct list_tag *next;
  point_t *p;
} list_t;

double *get_convex_minorant(count_t n, double *Gammas);

