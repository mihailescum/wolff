#ifndef WOLFF_GRAPHS_HIERARCHICAL_H
#define WOLFF_GRAPHS_HIERARCHICAL_H

#include "../wolff.hpp"

namespace wolff
{
    template <class vertex_prop = std::tuple<>, class edge_prop = double>
    class hierarchical_graph : public graph<vertex_prop, edge_prop>
    {
    public:
        unsigned R;

        hierarchical_graph() : graph<vertex_prop, edge_prop>(), R(0) {}

        /// @brief
        /// @param D Dimension of the graph
        /// @param L Size of box at highest recursion level
        /// @param R number of recursion levels
        hierarchical_graph(unsigned D, unsigned L, unsigned R) : graph<vertex_prop, edge_prop>(D, L, false), R(R)
        {
            if (R <= 0)
            {
                throw("At least one level of recursion is necessary");
            }
        }

        void init() override
        {
            // Need to declare `nv`, `ne` and `vertices`
            this->nv = get_depth_offset(R);
            this->ne = 0;
            this->vertices.resize(this->nv);

            unsigned nv_box = pow(this->L, this->D);
            for (unsigned depth = 0; depth < R; ++depth)
            {
                unsigned size_level = pow(this->L, depth * this->D);
                for (unsigned l = 0; l < size_level; ++l)
                {
                    unsigned offset = get_depth_offset(depth) + l * nv_box;
                    build_recursive_level(offset);
                    if (depth >= 1)
                    {
                        unsigned parent = get_depth_offset(depth - 1) + l;
                        connect_lower_level(parent, offset);
                    }

                    offset += nv_box;
                }
            }
            this->ne /= 2;
        }

        unsigned get_depth_offset(unsigned depth)
        {
            const unsigned nv_box = pow(this->L, this->D);
            return (pow(nv_box, depth + 1) - 1) / (nv_box - 1) - 1;
        }

        void connect_lower_level(unsigned parent, unsigned offset)
        {
            const unsigned nv_box = pow(this->L, this->D);

            for (unsigned i = 0; i < nv_box; ++i)
            {
                unsigned child = offset + i;
                edge_prop prop{1.0 / nv_box};
                add_edge(parent, child, prop);
                add_edge(child, parent, prop);
            }
        }

        void build_recursive_level(unsigned offset)
        {
            unsigned nv_level = pow(this->L, this->D);

            for (unsigned i = 0; i < nv_level; i++)
            {
                unsigned self = offset + i;
                this->vertices[self].ind = self;
                for (unsigned j = 0; j < this->D; j++)
                {
                    unsigned n1 = pow(this->L, j + 1) * (i / ((unsigned)pow(this->L, j + 1))) +
                                  fmod(i + pow(this->L, j), pow(this->L, j + 1));

                    if (this->touches_boundary(i, j) != 1) // We did not reach the end of the box
                    {
                        edge_prop prop{1.0};
                        add_edge(self, offset + n1, prop);
                    }

                    unsigned n2 = pow(this->L, j + 1) * (i / ((unsigned)pow(this->L, j + 1))) +
                                  fmod(pow(this->L, j + 1) + i - pow(this->L, j), pow(this->L, j + 1));

                    if (this->touches_boundary(i, j) != -1) // We did not reach the end of the box
                    {
                        edge_prop prop{1.0};
                        add_edge(self, offset + n2, prop);
                    }
                }
            }
        }
        void add_edge(unsigned self, unsigned neighbor, edge_prop prop)
        {
            typename graph<vertex_prop, edge_prop>::halfedge f(this->vertices[self], this->vertices[neighbor]);
            f.prop = prop;
            this->vertices[self].edges.push_back(f);
            this->ne++;
        }
    };
}

#endif