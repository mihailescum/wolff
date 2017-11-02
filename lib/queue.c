
#include "queue.h"

void stack_push(ll_t **q, uint32_t x) {
  ll_t *nq = malloc(sizeof(ll_t));
  nq->x = x;
  nq->next = *q;

  *q = nq;
}

void stack_push_d(dll_t **q, double x) {
  dll_t *nq = malloc(sizeof(dll_t));
  nq->x = x;
  nq->next = *q;

  *q = nq;
}

uint32_t stack_pop(ll_t **q) {
  ll_t *old_q = *q;

  *q = old_q->next;
  uint32_t x = old_q->x;

  free(old_q);

  return x;
}

double stack_pop_d(dll_t **q) {
  dll_t *old_q = *q;

  *q = old_q->next;
  double x = old_q->x;

  free(old_q);

  return x;
}

bool stack_contains(const ll_t *q, uint32_t x) {
  if (q == NULL) {
    return false;
  } else if (q->x == x) {
    return true;
  } else {
    return stack_contains(q->next, x);
  }
}
