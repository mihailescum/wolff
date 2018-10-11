
#pragma once

#include <inttypes.h>
#include <cmath>
#include <stdlib.h>
#include <vector>

#include "types.h"

class graph_t {
  public:
    v_t ne;
    v_t nv;
    std::vector<std::vector<v_t>> v_adj;
    std::vector<std::vector<double>> coordinate;

    graph_t(D_t D, L_t L);
    void add_ext();
};

