
#pragma once

#include <cmath>
#include <vector>

#include "types.h"

typedef enum lattice_t {
  SQUARE_LATTICE,
  DIAGONAL_LATTICE
} lattice_t;

class graph_t {
  public:
    v_t ne;
    v_t nv;
    std::vector<std::vector<v_t>> v_adj;
    std::vector<std::vector<double>> coordinate;

    graph_t(D_t D, L_t L, lattice_t lat = SQUARE_LATTICE);
    void add_ext();
};

