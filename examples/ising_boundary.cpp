
#include <getopt.h>
#include <iostream>
#include <chrono>
#include <vector>
#include <optional>

// #define WOLFF_USE_FINITE_STATES
#define WOLFF_SITE_DEPENDENCE
#include <wolff_models/ising.hpp>

#include "render2d_measurement.hpp"

using namespace wolff;

typedef graph<std::vector<std::optional<ising_t>>> box_graph;
typedef wolff::system<ising_t, ising_t, box_graph> sys;

void print(sys S)
{
    for (auto &v : S.G.vertices)
    {
        std::cout << v.ind << "\t";
        for (auto &e : v.edges)
        {
            std::cout << e.neighbor.ind << " ";
        }
        std::cout << "\n";
    }
}

int main(int argc, char *argv[])
{
    // set defaults
    unsigned N = (unsigned)1e4;
    unsigned D = 3;
    unsigned L = 16;
    double T = 2.0; // 2.26918531421;
    double H = 0.0;
    unsigned int window_size = 512;
    unsigned int frame_skip = 1;

    int opt;

    // take command line arguments
    while ((opt = getopt(argc, argv, "N:D:L:T:H:w:f:")) != -1)
    {
        switch (opt)
        {
        case 'N': // number of steps
            N = (unsigned)atof(optarg);
            break;
        case 'D': // dimension
            D = atoi(optarg);
            break;
        case 'L': // linear size
            L = atoi(optarg);
            break;
        case 'T': // temperature
            T = atof(optarg);
            break;
        case 'H': // external field
            H = atof(optarg);
            break;
        case 'w':
            window_size = atoi(optarg);
            break;
        case 'f':
            frame_skip = atoi(optarg);
            break;
        default:
            exit(EXIT_FAILURE);
        }
    }

    // define the spin-spin coupling
    std::function<double(const ising_t &, const ising_t &)> Z = [](const ising_t &s1, const ising_t &s2) -> double
    {
        return (double)(s1 * s2);
    };

    // define the spin-field coupling
    std::function<double(const box_graph::vertex &, const ising_t &)> B = [=](const box_graph::vertex &v, const ising_t &s) -> double
    {
        double result = H * s;
        for (auto &p : v.prop)
        {
            if (p.has_value())
            {
                result += Z(s, p.value());
            }
        }

        return result;
    };

    // initialize the lattice
    box_graph G(D, L, false);
    for (auto &v : G.vertices)
    {
        for (unsigned j = 0; j < D; j++)
        {
            int touches_boundary = G.touches_boundary(v.ind, j);
            if (touches_boundary == -1 || touches_boundary == 1)
            {
                const ising_t value = (G.coordinate_projection(v.ind, 0) < L / 4) ? false : true;
                v.prop.push_back(value);
            }
            else
            {
                v.prop.push_back({});
            }
        }
    }

    // initialize the system
    sys S(G, T, Z, B);
    // print(S);

    // initailze the measurement object
    std::function<Vector4f(const wolff::ising_t &, const sys &S)> color_conversion = [](const wolff::ising_t &s, const sys &S) -> Vector4f
    {
        Vector4f color;
        if (s.x == S.s0.x)
        {
            color = {1.0, 0.0, 0.0, 1.0};
        }
        else
        {
            color = {1.0, 1.0, 1.0, 1.0};
        }
        return color;
    };
    render2d_measurement<wolff::ising_t, wolff::ising_t, box_graph> A(0, 1, 0, L, color_conversion, window_size, "Ising", argc, argv);

    // initialize the random number generator
    auto seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
    std::mt19937 rng(seed);
    std::uniform_int_distribution<unsigned> dist(0, 1);

    // Random init
    for (unsigned i = 0; i < G.nv; i++)
    {
        S.s[i] = bool(dist(rng));
    }

    // run wolff N times
    S.run_wolff(N, gen_ising<box_graph>, A, rng);

    while (true)
    {
        A.post_cluster(-1, N, S);
    }

    // exit
    return 0;
}
