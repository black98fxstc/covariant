#include <cmath>
#include <array>
#include <vector>
#include <fftw3.h>

template <unsigned Dimension, unsigned Grid>
class Covariant
{
public:
    const unsigned dimension = Dimension;
    const unsigned grid = Grid;
    size_t stride[Dimension];
    const double delta = 1.0 / grid;
    typedef std::vector<float> Event;

private:
    static constexpr unsigned points = Grid + 1;
    static const size_t size = std::pow(points, Dimension);
    size_t _events = 0;
    float *_weight;
    float *_density;
    std::array<std::vector<float>, Dimension> _f;
    std::array<std::array<std::vector<float>, Dimension>, Dimension> _s;
    std::array<std::array<std::vector<float>, Dimension>, Dimension> _t;
    std::array<std::vector<float>, Dimension> _T;
    std::array<std::vector<float>, Dimension> _S;
    std::vector<float> _L;

public:
    struct Fiber
    {
    public:
        size_t id;
        size_t base;
        size_t stride;
        unsigned dimension;
    };

    bool
    event(Event &event)
    {
        size_t coordinate = 0;
        double rem[Dimension];
        for (unsigned i = 0; i < dimension; i++)
        {
            if (event[i] < 0.0f || event[i] >= 1.0f)
                return false;
            rem[i] = event[i] * grid;
            unsigned floor = static_cast<unsigned>(rem[i]);
            coordinate += floor * stride[i];
            rem[i] -= floor;
        }

        for (unsigned neighbor = 0; neighbor < (unsigned)(1 << dimension); neighbor++)
        {
            size_t offset = 0;
            double weight = 1.0;
            for (unsigned i = 0; i < Dimension; i++)
            {
                if (neighbor & (1 << i))
                {
                    offset += stride[i];
                    weight *= rem[i];
                }
                else
                {
                    weight *= 1.0f - rem[i];
                }
            }
            _weight[coordinate + offset] += weight;
        }
        ++_events;
        return true;
    }

    size_t events()
    {
        return _events;
    }

    void smooth(float percent)
    {
        if (percent <= 0.0f)
        _density = _weight;
    }

    void parameters()
    {
        for_each_fiber([this](Fiber& fiber) { this->basis_functions(fiber); });
        for_each_fiber([this](Fiber& fiber) { this->natural_parameters(fiber); });
        for (size_t i = 0; i < size; i++)
        {
            for (unsigned j = 0; j < Dimension; j++)
                for (unsigned k = 0; k < Dimension; k++)
                    _T[j][i] += _t[k][j][i];
            for (unsigned j = 0; j < Dimension; j++)
                _L[i] += _T[j][i];
        }
    }

    Covariant() {
        size_t s = 1;
        for (unsigned i = 0; i < dimension; i++)
        {
            stride[i] = s;
            s *= points;
        }
        _weight = (float*) fftw_malloc(sizeof(float) * size);
    }

private:
    void *basis_functions(Fiber &fiber)
    {
        double marginal = (_density[fiber.base] + _density[fiber.base + grid * fiber.stride]) / 2.0;
        for (size_t i = 1; i < grid; i++)
        {
            marginal += _density[fiber.base + i * fiber.stride];
        }
        for (size_t i = 0; i < points; i++)
        {
            _f[fiber.dimension][fiber.base + i * fiber.stride] = _density[fiber.base + i * fiber.stride] / marginal;
        }
        return nullptr;
    }

    void *natural_parameters(Fiber &fiber)
    {
        for (unsigned dim = 0; dim < Dimension; dim++)
        {
            float max = 0.0f;
            unsigned m = grid / 2;
            for (unsigned i = 0; i < points; i++)
            {
                if (_f[dim][fiber.base + i * fiber.stride] > max)
                {
                    max = _f[dim][fiber.base + i * fiber.stride];
                    m = i;
                }
            }
            double t;
            _s[dim][fiber.dimension][fiber.base + m * fiber.stride] = 0.0;
            for (unsigned i = m, j; i < points;)
            {
                j = i + 1;
                while (_f[dim][fiber.base + j * fiber.stride] <= 0.0 && j < grid)
                    j++;
                if (_t[dim][fiber.dimension][fiber.base + j * fiber.stride] <= 0.0)
                    t = 1.0 / delta / delta;
                else
                    t = -2.0 *
                        (std::log(_f[dim][fiber.base + j * fiber.stride]) - std::log(_f[dim][fiber.base + i * fiber.stride]) / delta / delta / (j - i) / (j - i) - _s[dim][fiber.dimension][fiber.base + i * fiber.stride] / delta / (j - i));
                while (i < j)
                {
                    _t[dim][fiber.dimension][fiber.base + i * fiber.stride] = t;
                    _s[dim][fiber.dimension][fiber.base + (i + 1) * fiber.stride] = -_t[dim][fiber.dimension][fiber.base + i * fiber.stride] * delta + _s[dim][fiber.dimension][fiber.base + i * fiber.stride];
                    i++;
                }
            }
            for (unsigned i = m, j; i > 0;)
            {
                j = i - 1;
                while (_f[dim][fiber.base + j * fiber.stride] == 0.0 && j > 0)
                    j--;
                if (_t[dim][fiber.dimension][fiber.base + j * fiber.stride] <= 0.0)
                    t = 1.0 / delta / delta;
                else
                    t = -2.0 *
                        (std::log(_f[dim][fiber.base + j * fiber.stride]) - std::log(_f[dim][fiber.base + i * fiber.stride]) / delta / delta / (j - i) / (j - i) - _s[dim][fiber.dimension][fiber.base + i * fiber.stride] / delta / (j - i));
                while (i > j)
                {
                    _t[dim][fiber.dimension][fiber.base + i * fiber.stride] = t;
                    _s[dim][fiber.dimension][fiber.base + (i + 1) * fiber.stride] = -_t[dim][fiber.dimension][fiber.base + i * fiber.stride] * delta + _s[dim][fiber.dimension][fiber.base + i * fiber.stride];
                    i--;
                }
            }
        }
        return nullptr;
    }

    void for_each_fiber(std::function<void(Fiber &)> func)
    {
        Fiber fiber;
        for (unsigned dim = 0; dim < Dimension; ++dim)
        {
            fiber.dimension = dim;
            fiber.stride = stride[dim];

            for (fiber.id = 0; fiber.id < size / points; fiber.id++)
            {
                size_t smaller = dim == 0 ? 1 : fiber.id % stride[dim];
                size_t larger = fiber.id / stride[dim];
                fiber.base = larger * points * smaller;
                func(fiber);
            }
        }
    }

};
