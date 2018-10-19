
#ifndef WOLFF_GRAPH_H
#define WOLFF_GRAPH_H

#include <cmath>
#include <vector>

namespace wolff {

#include "types.h"

class graph {
  public:
    D_t D;
    L_t L;
    v_t ne;
    v_t nv;
    std::vector<std::vector<v_t>> adjacency;
    std::vector<std::vector<double>> coordinate;

    graph();
    graph(D_t D, L_t L);

    void add_ghost();
};

}

#endif

