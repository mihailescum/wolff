
#include <getopt.h>
#include <iostream>
#include <chrono>
#include <optional>

#define WOLFF_SITE_DEPENDENCE
#include <wolff_models/height.hpp>
#include <wolff_models/dihedral_inf.hpp>
#include <wolff.hpp>

#include "overview_measurement.hpp"
#include "simple_measurement.hpp"
#include "render2d_measurement.hpp"

typedef graph<std::vector<std::optional<height_t<double>>>> box_graph;

int main(int argc, char *argv[])
{

    // set defaults
    unsigned N = (unsigned)250;
    unsigned D = 2;
    unsigned L = 1024;
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
    auto Z = [](const height_t<double> &s1, const height_t<double> &s2) -> double
    {
        return -pow(s1.x - s2.x, 2);
    };

    // define the spin-field coupling
    auto B = [&](const box_graph::vertex &v, const height_t<double> &s) -> double
    {
        double result = -H * pow(s.x, 2);

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
    double boundary_value = 0.0;
    for (auto &v : G.vertices)
    {
        for (unsigned j = 0; j < D; j++)
        {
            int touches_boundary = G.touches_boundary(v.ind, j);
            if (touches_boundary == -1 || touches_boundary == 1)
            {
                height_t<double> value(boundary_value);
                v.prop.push_back(value);
            }
            else
            {
                v.prop.push_back({});
            }
        }
    }

    // initialize the system
    wolff::system<dihedral_inf_t<double>, height_t<double>, box_graph> S(G, T, Z, B);

    bool odd_run = false;

    auto gen_R_IH = [&](std::mt19937 &r, const wolff::system<dihedral_inf_t<double>, height_t<double>, box_graph> &S, const box_graph::vertex &v) -> dihedral_inf_t<double>
    {
        dihedral_inf_t<double> rot;
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
            std::normal_distribution<double> dist(0.0, 1.0);
            rot.x = 2 * S.s[v.ind].x + dist(r);
        }

        odd_run = !odd_run;

        return rot;
    };

    // initailze the measurement object
    // overview_measurement<dihedral_inf_t<double>, height_t<double>, box_graph> A;
    simple_measurement<dihedral_inf_t<double>, height_t<double>, box_graph> A(S);

    // initialize the random number generator
    auto seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
    std::mt19937 rng(seed);

    // run wolff N times
    auto init_spins = [&](unsigned) -> height_t<double>
    {
        std::normal_distribution<double> initialization_distribution(0.0, 1.0);
        return height_t(initialization_distribution(rng) + boundary_value);
    };
    S.init_spins(init_spins);
    S.run_wolff(N, gen_R_IH, A, rng);

    // print the result of our measurements
    std::cout << "Wolff complete!\nThe average energy per site was " << A.avgE() / S.nv
              << ".\nThe average magnetization per site was " << A.avgM() / S.nv
              << ".\nThe average cluster size per site was " << A.avgC() / S.nv << ".\n";

    float spin_min = NAN;
    float spin_max = NAN;
    auto color_conversion = [&](unsigned index, const wolff::system<dihedral_inf_t<double>, height_t<double>, box_graph> &S) -> Vector4f
    {
        if (std::isnan(spin_min))
        {
            spin_min = S.s[0].x;
            spin_max = S.s[0].x;
            for (auto &s : S.s)
            {
                spin_min = std::min(static_cast<float>(s.x), spin_min);
                spin_max = std::max(static_cast<float>(s.x), spin_max);
            }
            std::cout << "Minimal spin: " << spin_min << "\t";
            std::cout << "Maximal spin: " << spin_max << "\n";
        }

        auto s = S.s[index];
        float r = (s.x - spin_min) / (spin_max - spin_min);
        if (r > 1.0)
            r = 1.0;
        else if (r < 0.0)
            r = 0.0;

        float g = 0.0;
        float b = 0.0;
        for (unsigned d = 0; d < D; d++)
        {
            int touches_boundary = G.touches_boundary(index, d);
            if (touches_boundary == 1)
            {
                g = 1.0;
            }
            else if (touches_boundary == -1)
            {
                b = 1.0;
            }
        }

        Vector4f color{r, g, b, 1.0};
        return color;
    };
    unsigned int window_size = std::max(512U, L);
    render2d_measurement<dihedral_inf_t<double>, height_t<double>, box_graph> render(0, 1, 0, L, color_conversion, window_size, "Continuous Gaussian", argc, argv);

    render.show(S);

    // exit
    return 0;
}
