#include <cmath>
#include <array>
#include <vector>
#include <fftw3.h>
#include <assert.h>

template <unsigned Dimension, unsigned Grid>
class Covariant
{
public:
    const unsigned dimension = Dimension;
    const unsigned grid = Grid;
    size_t stride[Dimension];
    int points[Dimension];
    const double delta = 1.0 / grid;
    typedef std::vector<float> Event;
    struct Fiber
    {
    public:
        size_t id;
        size_t base;
        size_t stride;
        unsigned d;
    };

private:
    static const size_t size = std::pow(Grid + 1, Dimension);
    size_t _events = 0;
    float *_weight;
    float *_density;
    fftw_r2r_kind kind[Dimension];
    void *DCT = nullptr;
    std::array<std::array<float, size>, Dimension> _f;
    std::array<std::array<std::array<float, size>, Dimension>, Dimension> _s;
    std::array<std::array<std::array<float, size>, Dimension>, Dimension> _t;
    std::array<std::array<float, size>, Dimension> _T;
    std::array<std::array<float, size>, Dimension> _S;
    std::array<float, size> _L;

public:
    bool
    event(Event &event)
    {
        size_t x = 0;
        double rem[Dimension];
        for (unsigned i = 0; i < dimension; i++)
        {
            if (event[i] < 0.0f || event[i] >= 1.0f)
                return false;
            rem[i] = event[i] * grid;
            unsigned floor = static_cast<unsigned>(rem[i]);
            x += floor * stride[i];
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
            _weight[x + offset] += weight;
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
        for_each_fiber([this](Fiber &fiber)
                       { this->basis_functions(fiber); });
        for_each_fiber([this](Fiber &fiber)
                       { this->natural_parameters(fiber); });
        for (size_t x = 0; x < size; x++)
        {
            for (unsigned j = 0; j < Dimension; j++)
                for (unsigned i = 0; i < Dimension; i++)
                    _T[j][x] += _t[i][j][x];
            for (unsigned j = 0; j < Dimension; j++)
                _L[x] += _T[j][x];
        }
    }

    Covariant()
    {
        size_t s = 1;
        for (unsigned i = 0; i < dimension; i++)
        {
            points[i] = Grid + 1;
            stride[i] = s;
            kind[i] = FFTW_REDFT00;
            s *= points[i];
        }
        _weight = (float *)fftw_malloc(sizeof(float) * size);
        _density = (float *)fftw_malloc(sizeof(float) * size);
        DCT = (void *)fftwf_plan_r2r(Dimension, points, _weight, _density, kind, 0);
        assert(DCT);
    }

    ~Covariant()
    {
        fftwf_destroy_plan((fftwf_plan)DCT);
    };

private:
    void *basis_functions(Fiber &fiber)
    {
        double marginal = (_density[fiber.base] + _density[fiber.base + grid * fiber.stride]) / 2.0;
        for (size_t i = 1; i < grid; i++)
        {
            marginal += _density[fiber.base + i * fiber.stride];
        }
        for (int i = 0; i < points[fiber.d]; i++)
        {
            _f[fiber.d][fiber.base + i * fiber.stride] = _density[fiber.base + i * fiber.stride] / marginal;
        }
        return nullptr;
    }

    void *natural_parameters(Fiber &fiber)
    {
        for (unsigned i = 0; i < Dimension; i++)
        {
            float max = 0.0f;
            unsigned m = grid / 2;
            for (int j = 0; j < points[fiber.d]; j++)
            {
                if (_f[i][fiber.base + j * fiber.stride] > max)
                {
                    max = _f[i][fiber.base + j * fiber.stride];
                    m = j;
                }
            }
            double t;
            _s[i][fiber.d][fiber.base + m * fiber.stride] = 0.0;
            for (int k = m, j; k < points[fiber.d];)
            {
                j = k + 1;
                while (_f[i][fiber.base + j * fiber.stride] <= 0.0 && j < points[fiber.d] - 1)
                    j++;
                if (_t[i][fiber.d][fiber.base + j * fiber.stride] <= 0.0)
                    t = 1.0 / delta / delta;
                else
                    t = -2.0 * (std::log(_f[i][fiber.base + j * fiber.stride]) - std::log(_f[i][fiber.base + k * fiber.stride]) / delta / delta / (j - k) / (j - k) 
                    - _s[i][fiber.d][fiber.base + k * fiber.stride] / delta / (j - k));
                while (k < j)
                {
                    _t[i][fiber.d][fiber.base + k * fiber.stride] = t;
                    _s[i][fiber.d][fiber.base + (k + 1) * fiber.stride] = -_t[i][fiber.d][fiber.base + k * fiber.stride] * delta + _s[i][fiber.d][fiber.base + k * fiber.stride];
                    k++;
                }
            }
            for (unsigned k = m, j; k > 0;)
            {
                j = k - 1;
                while (_f[i][fiber.base + j * fiber.stride] == 0.0 && j > 0)
                    j--;
                if (_t[i][fiber.d][fiber.base + j * fiber.stride] <= 0.0)
                    t = 1.0 / delta / delta;
                else
                    t = -2.0 * (std::log(_f[i][fiber.base + j * fiber.stride]) - std::log(_f[i][fiber.base + k * fiber.stride]) / delta / delta / (j - k) / (j - k) 
                    - _s[i][fiber.d][fiber.base + k * fiber.stride] / delta / (j - k));
                while (k > j)
                {
                    _t[i][fiber.d][fiber.base + k * fiber.stride] = t;
                    _s[i][fiber.d][fiber.base + (k + 1) * fiber.stride] = -_t[i][fiber.d][fiber.base + k * fiber.stride] * delta + _s[i][fiber.d][fiber.base + k * fiber.stride];
                    k--;
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
            fiber.d = dim;
            fiber.stride = stride[dim];

            for (fiber.id = 0; fiber.id < size / points[fiber.d]; fiber.id++)
            {
                size_t smaller = dim == 0 ? 1 : fiber.id % stride[dim];
                size_t larger = fiber.id / stride[dim];
                fiber.base = larger * points[fiber.d] * smaller;
                func(fiber);
            }
        }
    }
};
