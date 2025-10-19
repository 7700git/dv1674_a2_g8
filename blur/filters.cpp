/*
Author: David Holmqvist <daae19@student.bth.se>
*/

#include "filters.hpp"
#include "matrix.hpp"
#include "ppm.hpp"
#include <cmath>
#include <cstdef>

namespace Filter
{

    namespace Gauss
    {
        void get_weights(int n, double *weights_out)
        {
            for (auto i{0}; i <= n; i++)
            {
                double x{static_cast<double>(i) * max_x / n};
                weights_out[i] = exp(-x * x * pi);
            }
        }
    }

    Matrix blur(Matrix m, const int radius)
    {
        Matrix scratch{PPM::max_dimension};
        auto dst{m};
	//move the weight calc to outside of nested loops to exponentially decrease calls to this function
        double w[Gauss::max_radius]{};
        Gauss::get_weights(radius, w);

	std::size_t width = dst.get_x_size();
	std::size_t height = dst.get_y_size();
	//the loops are based on scanning columns first which uses more cache, hence flipping the loop inside out is more efficient
        for (std::size_t y=0; y < height; ++y)
        {
            for (std::size_t x=0; x < width; ++x)
            {

                // unsigned char Matrix::r(unsigned x, unsigned y) const
                // {
                //     return R[y * x_size + x];
                // }

                auto r{w[0] * dst.r(x, y)}, g{w[0] * dst.g(x, y)}, b{w[0] * dst.b(x, y)}, n{w[0]};

                for (int wi=1; wi <= radius; ++wi)
                {
                    double wc = w[wi];
                    int x2 = static_cast<int>(x) - wi;
                    if (x2 >= 0)
                    {
                        r += wc * dst.r(x2, y);
                        g += wc * dst.g(x2, y);
                        b += wc * dst.b(x2, y);
                        n += wc;
                    }
                    x2 = static_cast<int>x + wi;
                    if (x2 < static_cast<int>(width))
                    {
                        r += wc * dst.r(x2, y);
                        g += wc * dst.g(x2, y);
                        b += wc * dst.b(x2, y);
                        n += wc;
                    }
                }
                scratch.r(x, y) = r / n;
                scratch.g(x, y) = g / n;
                scratch.b(x, y) = b / n;
            }
        }
	//same alteration in these loops as above
        //double w[Gauss::max_radius]{};
        //Gauss::get_weights(radius, w);
	//same change as above, flipping the loop matrix inside out
        for (std::size_t y = 0; y < height; ++y)
        {
            for (std::size_t x = 0; x < width; ++x)
            {

                auto r{w[0] * scratch.r(x, y)}, g{w[0] * scratch.g(x, y)}, b{w[0] * scratch.b(x, y)}, n{w[0]};

                for (int wi = 1; wi <= radius; ++wi)
                {
                    double wc = w[wi];
                    int y2 = static_cast<int>(y) - wi;
                    if (y2 >= 0)
                    {
                        r += wc * scratch.r(x, y2);
                        g += wc * scratch.g(x, y2);
                        b += wc * scratch.b(x, y2);
                        n += wc;
                    }
                    y2 = statis_cast<int>(y) + wi;
                    if (y2 < static_cast<int>(height))
                    {
                        r += wc * scratch.r(x, y2);
                        g += wc * scratch.g(x, y2);
                        b += wc * scratch.b(x, y2);
                        n += wc;
                    }
                }
                dst.r(x, y) = r / n;
                dst.g(x, y) = g / n;
                dst.b(x, y) = b / n;
            }
        }

        return dst;
    }

}
