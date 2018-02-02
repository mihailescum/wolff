
#include "tree.h"

node_t *tree_createNode(v_t value, v_t level, node_t *left, node_t *right) {
  node_t *node = (node_t *)malloc(sizeof(node_t));

  node->value = value;
  node->level = level;
  node->left = left;
  node->right = right;

  return node;
}

void tree_freeNode(node_t *node) {
  if (node != NULL) {
    tree_freeNode(node->left);
    tree_freeNode(node->right);
    free(node);
  }
}

void tree_skew(node_t **T) {
  if (*T != NULL) {
    if ((*T)->left != NULL) {
      if ((*T)->left->level == (*T)->level) {
        node_t *L = (*T)->left;
        (*T)->left = L->right;
        L->right = (*T);

        *T = L;
      }
    }
  }
}

void tree_split(node_t **T) {
  if (*T != NULL) {
    if ((*T)->right != NULL) {
      if ((*T)->right->right != NULL) {
        if ((*T)->level == (*T)->right->right->level) {
          node_t *R = (*T)->right;
          (*T)->right = R->left;
          R->left = *T;
          R->level = R->level + 1;
          *T = R;
        }
      }
    }
  }
}

void tree_insert(node_t **T, v_t x) {
  if (*T == NULL) {
    *T = tree_createNode(x, 1, NULL, NULL);
  } else if (x < (*T)->value) {
    node_t *L = (*T)->left;
    tree_insert(&L, x);
    (*T)->left = L;
  } else if (x > (*T)->value) {
    node_t *R = (*T)->right;
    tree_insert(&R, x);
    (*T)->right = R;
  }

  tree_skew(T);
  tree_split(T);
}

bool tree_contains(node_t *T, v_t x) {
  if (T == NULL) {
    return false;
  } else if (x < T->value) {
    return tree_contains(T->left, x);
  } else if (x > T->value) {
    return tree_contains(T->right, x);
  } else {
    return true;
  }
}

