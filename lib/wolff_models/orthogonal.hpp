
#ifndef WOLFF_MODELS_ORTHOGONAL_H
#define WOLFF_MODELS_ORTHOGONAL_H

#include <random>
#include <cmath>

#include <wolff.hpp>
#include "vector.hpp"

namespace wolff
{

    template <unsigned q, class T>
    class orthogonal_t : public std::array<std::array<T, q>, q>
    {
    public:
        bool is_reflection;

        orthogonal_t() : is_reflection(false)
        {
            for (unsigned i = 0; i < q; i++)
            {
                (*this)[i].fill(0);
                (*this)[i][i] = (T)1;
            }
        }

        vector_t<q, T> act(const vector_t<q, T> &v) const
        {
            vector_t<q, T> v_rot;
            v_rot.fill(0);

            if (is_reflection)
            {
                double prod = 0;
                for (unsigned i = 0; i < q; i++)
                {
                    prod += v[i] * (*this)[0][i];
                }
                for (unsigned i = 0; i < q; i++)
                {
                    v_rot[i] = v[i] - 2 * prod * (*this)[0][i];
                }
            }
            else
            {
                for (unsigned i = 0; i < q; i++)
                {
                    for (unsigned j = 0; j < q; j++)
                    {
                        v_rot[i] += (*this)[i][j] * v[j];
                    }
                }
            }

            return v_rot;
        }

        orthogonal_t<q, T> act(const orthogonal_t<q, T> &m) const
        {
            orthogonal_t<q, T> m_rot;

            m_rot.is_reflection = false;

            if (is_reflection)
            {
                for (unsigned i = 0; i < q; i++)
                {
                    double akOki = 0;

                    for (unsigned k = 0; k < q; k++)
                    {
                        akOki += (*this)[0][k] * m[k][i];
                    }

                    for (unsigned j = 0; j < q; j++)
                    {
                        m_rot[j][i] = m[j][i] - 2 * akOki * (*this)[0][j];
                    }
                }
            }
            else
            {
                for (unsigned i = 0; i < q; i++)
                {
                    m_rot[i].fill(0);
                    for (unsigned j = 0; j < q; j++)
                    {
                        for (unsigned k = 0; k < q; k++)
                        {
                            m_rot[i][j] += (*this)[i][j] * m[j][k];
                        }
                    }
                }
            }

            return m_rot;
        }

        vector_t<q, T> act_inverse(const vector_t<q, T> &v) const
        {
            if (is_reflection)
            {
                return this->act(v); // reflections are their own inverse
            }
            else
            {
                vector_t<q, T> v_rot;
                v_rot.fill(0);

                for (unsigned i = 0; i < q; i++)
                {
                    for (unsigned j = 0; j < q; j++)
                    {
                        v_rot[i] += (*this)[j][i] * v[j];
                    }
                }

                return v_rot;
            }
        }

        vector_t<q, T> act_inverse(const orthogonal_t<q, T> &m) const
        {
            if (is_reflection)
            {
                return this->act(m); // reflections are their own inverse
            }
            else
            {
                orthogonal_t<q, T> m_rot;
                m_rot.is_reflection = false;

                for (unsigned i = 0; i < q; i++)
                {
                    m_rot[i].fill(0);
                    for (unsigned j = 0; j < q; j++)
                    {
                        for (unsigned k = 0; k < q; k++)
                        {
                            m_rot[i][j] += (*this)[j][i] * m[j][k];
                        }
                    }
                }

                return m_rot;
            }
        }
    };

    template <unsigned q, class G_t>
    orthogonal_t<q, double> generate_rotation_uniform(std::mt19937 &r, const system<orthogonal_t<q, double>, vector_t<q, double>, G_t> &, const typename G_t::vertex &)
    {
        std::normal_distribution<double> dist(0.0, 1.0);
        orthogonal_t<q, double> ptr;
        ptr.is_reflection = true;

        double v2 = 0;

        for (unsigned i = 0; i < q; i++)
        {
            ptr[0][i] = dist(r);
            v2 += ptr[0][i] * ptr[0][i];
        }

        double mag_v = sqrt(v2);

        for (unsigned i = 0; i < q; i++)
        {
            ptr[0][i] /= mag_v;
        }

        return ptr;
    }

    template <unsigned q, class G_t>
    orthogonal_t<q, double> generate_rotation_perturbation(std::mt19937 &r, const system<orthogonal_t<q, double>, vector_t<q, double>, G_t> &S, const typename G_t::vertex &v0, double epsilon, unsigned int n)
    {
        std::normal_distribution<double> dist(0.0, 1.0);
        orthogonal_t<q, double> m;
        m.is_reflection = true;

        vector_t<q, double> v;

        if (n > 1)
        {
            std::uniform_int_distribution<unsigned int> udist(0, n);
            unsigned int rotation = udist(r);

            double cosr = cos(2 * M_PI * rotation / (double)n / 2.0);
            double sinr = sin(2 * M_PI * rotation / (double)n / 2.0);

            v[0] = S.s[v0.ind][0] * cosr - S.s[v0.ind][1] * sinr;
            v[1] = S.s[v0.ind][1] * cosr + S.s[v0.ind][0] * sinr;

            for (unsigned i = 2; i < q; i++)
            {
                v[i] = S.s[v0.ind][i];
            }
        }
        else
        {
            v = S.s[v0.ind];
        }

        double m_dot_v = 0;

        for (unsigned i = 0; i < q; i++)
        {
            m[0][i] = dist(r); // create a random vector
            m_dot_v += m[0][i] * v[i];
        }

        double v2 = 0;

        for (unsigned i = 0; i < q; i++)
        {
            m[0][i] = m[0][i] - m_dot_v * v[i]; // find the component orthogonal to v
            v2 += pow(m[0][i], 2);
        }

        double mag_v = sqrt(v2);

        for (unsigned i = 0; i < q; i++)
        {
            m[0][i] /= mag_v; // normalize
        }

        v2 = 0;

        double factor = epsilon * dist(r);

        for (unsigned i = 0; i < q; i++)
        {
            m[0][i] += factor * v[i]; // perturb orthogonal vector in original direction
            v2 += pow(m[0][i], 2);
        }

        mag_v = sqrt(v2);

        for (unsigned i = 0; i < q; i++)
        {
            m[0][i] /= mag_v; // normalize
        }

        return m;
    }

}

#endif
