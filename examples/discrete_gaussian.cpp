
#include <getopt.h>
#include <iostream>
#include <chrono>
#include <optional>

#define WOLFF_SITE_DEPENDENCE
#include <wolff_models/height.hpp>
#include <wolff_models/dihedral_inf.hpp>
#include <wolff.hpp>

#include "simple_measurement.hpp"
#include "render2d_measurement.hpp"

typedef graph<std::vector<std::optional<height_t<long long>>>> box_graph;

int main(int argc, char *argv[])
{

    // set defaults
    unsigned N = (unsigned)1e4;
    unsigned D = 2;
    unsigned L = 128;
    double T = 0.8;
    double H = 0.0;

    int opt;

    // take command line arguments
    while ((opt = getopt(argc, argv, "N:D:L:T:H:")) != -1)
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
        default:
            exit(EXIT_FAILURE);
        }
    }

    // define the spin-spin coupling
    std::function<double(const height_t<long long> &, const height_t<long long> &)> Z = [](const height_t<long long> &s1, const height_t<long long> &s2) -> double
    {
        return -pow(s1.x - s2.x, 2);
    };

    // define the spin-field coupling
    auto B = [&](const box_graph::vertex &v, const height_t<long long> &s) -> double
    {
        double result = -H * pow(s.x - 10, 2);

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
    long long boundary_value = 0;
    for (auto &v : G.vertices)
    {
        for (unsigned j = 0; j < D; j++)
        {
            int touches_boundary = G.touches_boundary(v.ind, j);
            if (touches_boundary == -1 || touches_boundary == 1)
            {
                height_t<long long> value(boundary_value);
                v.prop.push_back(value);
            }
            else
            {
                v.prop.push_back({});
            }
        }
    }

    // initialize the system
    wolff::system<dihedral_inf_t<long long>, height_t<long long>, box_graph> S(G, T, Z, B);

    bool odd_run = false;

    auto gen_R_IH = [&](std::mt19937 &r, const wolff::system<dihedral_inf_t<long long>, height_t<long long>, box_graph> &S, const box_graph::vertex &v) -> dihedral_inf_t<long long>
    {
        dihedral_inf_t<long long> rot;
        rot.is_reflection = true;

        if (odd_run)
        {
            std::uniform_int_distribution<unsigned> dist(0, S.nv - 2);
            unsigned j = v.ind;

            // while (S.s[j].x == S.s[i0].x) {
            unsigned tmp = dist(r);

            if (tmp < v.ind)
            {
                j = tmp;
            }
            else
            {
                j = tmp + 1;
            }
            //}

            rot.x = 2 * S.s[j].x;
        }
        else
        {
            std::uniform_int_distribution<int> dist(0, 1);
            int j = dist(r);
            if (j)
            {
                rot.x = 2 * S.s[v.ind].x + 1;
            }
            else
            {
                rot.x = 2 * S.s[v.ind].x - 1;
            }
        }

        odd_run = !odd_run;

        return rot;
    };

    // initailze the measurement object
    simple_measurement A(S);

    // initialize the random number generator
    auto seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
    std::mt19937 rng(seed);

    // run wolff N times
    S.run_wolff(N, gen_R_IH, A, rng);

    // print the result of our measurements
    std::cout << "Wolff complete!\nThe average energy per site was " << A.avgE() / S.nv
              << ".\nThe average magnetization per site was " << A.avgM() / S.nv
              << ".\nThe average cluster size per site was " << A.avgC() / S.nv << ".\n";

    bool spin_set = false;
    long long spin_min = 0;
    long long spin_max = 0;
    auto color_conversion = [&](const height_t<long long> &s, const wolff::system<dihedral_inf_t<long long>, height_t<long long>, box_graph> &S) -> Vector4f
    {
        if (!spin_set)
        {
            spin_set = true;

            spin_min = S.s[0].x;
            spin_max = S.s[0].x;
            for (auto &s : S.s)
            {
                spin_min = std::min(static_cast<long long>(s.x), spin_min);
                spin_max = std::max(static_cast<long long>(s.x), spin_max);
            }
            std::cout << "Minimal spin: " << spin_min << "\t";
            std::cout << "Maximal spin: " << spin_max << "\n";
        }

        float r = float(s.x - spin_min) / (spin_max - spin_min);
        if (r > 1.0)
            r = 1.0;
        else if (r < 0.0)
            r = 0.0;

        float g = r; // cos(s * sqrt(T));
        float b = r;
        Vector4f color{r, g, b, 1.0};
        return color;
    };
    unsigned int window_size = std::max(512U, L);
    render2d_measurement<dihedral_inf_t<long long>, height_t<long long>, box_graph> render(0, 1, 0, L, color_conversion, window_size, "Discrete Gaussian", argc, argv);

    render.show(S);

    // exit
    return 0;
}
