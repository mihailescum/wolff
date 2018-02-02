
#pragma once

#include <inttypes.h>
#include <stdlib.h>
#include <stdbool.h>

#include "types.h"

typedef struct node_t {
  v_t value;
  v_t level;
  struct node_t *left;
  struct node_t *right;
} node_t;

void tree_insert(node_t **T, v_t x);

bool tree_contains(node_t *T, v_t x);

void tree_freeNode(node_t *T);

