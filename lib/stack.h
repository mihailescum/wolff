
#pragma once

#include <inttypes.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

#include "types.h"


#ifdef __cplusplus
extern "C" {
#endif

typedef struct ll_tag {
  v_t x;
  struct ll_tag *next;
} ll_t;

typedef struct dll_tag {
  double x;
  struct dll_tag *next;
} dll_t;

void stack_push(ll_t **q, v_t x);
void stack_push_d(dll_t **q, double x);

v_t stack_pop(ll_t **q);
double stack_pop_d(dll_t **q);

#ifdef __cplusplus
}
#endif

