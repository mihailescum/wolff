
#include "convex.h"

double slope(point_t *P, point_t *Q) {
  return (Q->y - P->y) / ((double)(Q->x) - (double)(P->x));
}

double *get_convex_minorant(count_t n, double *Gammas) {
  if (n < 2) {
    return Gammas;
  }

  list_t *L = (list_t *)calloc(1, sizeof(list_t));
  L->p = (point_t *)calloc(1, sizeof(point_t));
  L->p->x = 0;
  L->p->y = Gammas[0];

  list_t *pos = L;

  for (count_t i = 1; i < n; i++) {
    pos->next = (list_t *)calloc(1, sizeof(list_t));
    pos->next->p = (point_t *)calloc(1, sizeof(point_t));
    pos->next->p->x = i;
    pos->next->p->y = Gammas[i];
    pos->next->prev = pos;
    pos = pos->next;
  }

  pos->next = (list_t *)calloc(1, sizeof(list_t));
  pos->next->p = (point_t *)calloc(1, sizeof(point_t));
  pos->next->p->x = n;
  pos->next->p->y = 0;
  pos->next->prev = pos;

  list_t *X = L;
  list_t *Y = L->next;
  list_t *Z = Y->next;

  while (true) {
    if (slope(X->p, Y->p) <= slope(Y->p, Z->p)) {
      X = Y;
      Y = Z;
      if (Z->next == NULL) {
        break;
      } else {
        Z = Z->next;
      }
    } else {
      Y->prev->next = Y->next;
      Y->next->prev = Y->prev;
      free(Y->p);
      free(Y);
      if (X->prev != NULL) {
        Y = X;
        X = X->prev;
      } else {
        if (Z->next != NULL) {
          Y = Z;
          Z = Z->next;
        } else {
          break;
        }
      }
    }
  }

  pos = L;

  double *g = (double *)calloc(n + 1, sizeof(double));
  double rho = 0;

  for (count_t i = 0; i < n + 1; i++) {
    if (i > pos->next->p->x) {
      pos = pos->next;
    }

    g[i] = pos->p->y + ((double)i - (double)(pos->p->x)) *  (pos->next->p->y - pos->p->y) / ((double)(pos->next->p->x) - (double)(pos->p->x));

    if (i <n) {
      if (Gammas[i] - g[i] > rho) {
        rho = Gammas[i] - g[i];
      }
    } else {
      if (0 - g[i] > rho) {
        rho = 0 - g[i];
      }
    }
  }

  for (count_t i = 0; i < n + 1; i++) {
    g[i] += rho / 2;
  }

  while (L != NULL) {
    free(L->p);
    list_t *L_save = L;
    L = L->next;
    free(L);
  }

  return g;
}

