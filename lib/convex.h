
#pragma once

#include <inttypes.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
  uint64_t x;
  double y;
} point_t;

typedef struct list_tag {
  struct list_tag *prev;
  struct list_tag *next;
  point_t *p;
} list_t;

double *get_convex_minorant(uint64_t n, double *Gammas);

