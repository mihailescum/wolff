#ifndef __WOLFF___EXPORT_MEASUREMENT_HPP__
#define __WOLFF___EXPORT_MEASUREMENT_HPP__

#include <string>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <fstream>

#include <wolff.hpp>

template <class R_t, class X_t, class G_t>
class export_measurement : public wolff::measurement<R_t, X_t, G_t>
{
private:
    const wolff::system<R_t, X_t, G_t> &S;
    int n;

public:
    export_measurement(wolff::system<R_t, X_t, G_t> &S) : S(S), n(0) {}

    void pre_cluster(unsigned, unsigned, const wolff::system<R_t, X_t, G_t> &, const typename G_t::vertex &, const R_t &) override
    {
    }

    void plain_bond_visited(const wolff::system<R_t, X_t, G_t> &, const typename G_t::halfedge &, const X_t &, double dE) override
    {
    }

    void plain_site_transformed(const wolff::system<R_t, X_t, G_t> &S, const typename G_t::vertex &v, const X_t &si_new) override
    {
    }

    void post_cluster(unsigned n, unsigned N, const wolff::system<R_t, X_t, G_t> &S) override
    {
        int freq = (N > 1000) ? N / 1000 : 1;
        if (n <= N && n % freq == 0)
        {
            std::cout << "Iteration: " << n + 1 << "/" << N << "\n";
        }
        this->n++;
    }

    void export_results(const std::string &filename)
    {
        std::stringstream result_array;
        result_array << std::setprecision(4);
        for (const auto &s : S.s)
        {
            for (const auto &e : s)
            {
                result_array << e << " ";
            }
            // result_array << "\n";
        }

        std::ofstream outFile(filename);
        outFile << result_array.rdbuf();
        outFile.close();
    }
};

#endif