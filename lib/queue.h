
#pragma once

#include <inttypes.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

typedef struct ll_tag {
  uint32_t x;
  struct ll_tag *next;
} ll_t;

void stack_push(ll_t **q, uint32_t x);

uint32_t stack_pop(ll_t **q);

bool stack_contains(const ll_t *q, uint32_t x);
