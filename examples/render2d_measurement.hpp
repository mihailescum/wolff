#ifndef __WOLFF___RENDER2D_MEASUREMENT_HPP__
#define __WOLFF___RENDER2D_MEASUREMENT_HPP__

#include <iostream>

#include "progress_measurement.hpp"

#include "rendering.hpp"

using namespace wolff;

template <class R_t, class X_t, class G_t>
class render2d_measurement : public progress_measurement<R_t, X_t, G_t>
{
private:
    unsigned d1;
    unsigned d2;
    unsigned pivot;
    unsigned L;

    Renderer renderer;
    std::function<Vector4f(unsigned, const wolff::system<R_t, X_t, G_t> &)> color_conversion;

public:
    render2d_measurement(unsigned d1, unsigned d2, unsigned pivot, unsigned int L, std::function<Vector4f(unsigned, const wolff::system<R_t, X_t, G_t> &)> color_conversion, unsigned int window_size, const std::string title, int argc, char *argv[])
        : progress_measurement<R_t, X_t, G_t>(), d1(d1), d2(d2), pivot(pivot), L(L), color_conversion(color_conversion)
    {
        renderer.create(window_size, title);
        renderer.set_projection({0.0, 0.0, static_cast<float>(L), static_cast<float>(L)});

        glClearColor(0.0, 0.0, 0.0, 0.0);
    }

    void post_cluster(unsigned n, unsigned N, const wolff::system<R_t, X_t, G_t> &S) override
    {
        progress_measurement<R_t, X_t, G_t>::post_cluster(n, N, S);

        if (!renderer.is_running())
        {
            exit(EXIT_SUCCESS);
        }

        glClear(GL_COLOR_BUFFER_BIT);
        for (unsigned x = 0; x < L; x++)
        {
            for (unsigned y = 0; y < L; y++)
            {
                unsigned i = pivot + x * pow(L, d1) + y * pow(L, d2);

                Vector4f color = color_conversion(i, S);
                Vector2f position{static_cast<float>(x), static_cast<float>(y)};
                Vector2f scale{1.0f, 1.0f};
                renderer.draw_rect(position, scale, color);
            }
        }
        renderer.end_draw();
    }

    void show(const wolff::system<R_t, X_t, G_t> &S)
    {
        post_cluster(-1, 0, S);
        while (renderer.is_running())
        {
            glfwWaitEvents();
        }
    }
};

#endif