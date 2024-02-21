#ifndef __WOLFF___PROGRESS_MEASUREMENT_HPP__
#define __WOLFF___PROGRESS_MEASUREMENT_HPP__

#include <iostream>

#include "wolff.hpp"

using namespace wolff;

template <class R_t, class X_t, class G_t>
class progress_measurement : public measurement<R_t, X_t, G_t>
{

public:
    progress_measurement() : measurement<R_t, X_t, G_t>() {}

    void post_cluster(unsigned n, unsigned N, const wolff::system<R_t, X_t, G_t> &) override
    {
        if (n <= N)
        {
            std::cout << "Iteration: " << n + 1 << "/" << N << "\n";
        }
    }
};

#endif