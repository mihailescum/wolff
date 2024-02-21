#ifndef __WOLFF___OVERVIEW_MEASUREMENT_HPP__
#define __WOLFF___OVERVIEW_MEASUREMENT_HPP__

#include <wolff.hpp>

using namespace wolff;

template <class R_t, class X_t, class G_t>
class overview_measurement : public measurement<R_t, X_t, G_t>
{
public:
    overview_measurement() : measurement<R_t, X_t, G_t>()
    {
    }

    void post_cluster(unsigned n, unsigned N, const wolff::system<R_t, X_t, G_t> &S) override
    {
        double spin_average = 0.0;
        double spin_variance = 0.0;
        for (auto &s : S.s)
        {
            spin_average += s.x;
            spin_variance += s.x * s.x;
        }
        spin_average /= S.nv;
        spin_variance = spin_variance / S.nv - spin_average * spin_average;

        std::cout << "Iteration: " << n << "/" << N << "\t";
        std::cout << "Average spin: " << spin_average << "\t";
        std::cout << "Spin std: " << sqrt(spin_variance) << "\n";
    }
};

#endif