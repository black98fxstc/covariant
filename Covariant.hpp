#include <cmath>
#include <array>
#include <vector>
#include <fftw3.h>
#include <assert.h>

template <unsigned Dimension>
class Covariant
{
public:
    typedef std::vector<float> Event;
    const unsigned dimension = Dimension;

private:
    struct Fiber
    {
        size_t id;
        size_t base;
        size_t stride;
        unsigned d;
    };

    size_t _size;
    size_t _events = 0;
    float *_weight;
    float *_density;
    size_t stride[Dimension];
    int points[Dimension];
    fftw_r2r_kind kind[Dimension];
    void *DCT = nullptr;
    std::array<std::vector<float>, Dimension> _f;
    std::array<std::array<std::vector<float>, Dimension>, Dimension> _s;
    std::array<std::array<std::vector<float>, Dimension>, Dimension> _t;
    std::array<std::vector<float>, Dimension> _T;
    std::array<std::vector<float>, Dimension> _S;
    std::vector<float> _L;
    std::vector<float> _P;

public:
    size_t size() const
    {
        return _size;
    }

    size_t events()
    {
        return _events;
    }

    const float& L(size_t x) const
    {
        return _L[x];
    }

    const float& P(size_t x) const
    {
        return _P[x];
    }

    const float& f(unsigned i, size_t x) const
    {
        return _f[i][x];
    }

    const float& s(unsigned i, unsigned j, size_t x) const
    {
        return _s[i][j][x];
    }

    const float& t(unsigned i, unsigned j, size_t x) const
    {
        return _t[i][j][x];
    }

    bool event(Event &event)
    {
        size_t x = 0;
        double rem[Dimension];
        for (unsigned i = 0; i < Dimension; i++)
        {
            if (event[i] < 0.0f || event[i] >= 1.0f)
                return false;
            rem[i] = event[i] * (points[i] - 1);
            unsigned floor = static_cast<unsigned>(rem[i]);
            x += floor * stride[i];
            rem[i] -= floor;
        }

        for (unsigned neighbor = 0; neighbor < (unsigned)(1 << Dimension); neighbor++)
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

    void parameters(float percent = 1.0f)
    {
        if (percent > 0.0f)
        {
            float *cosine = (float *)fftw_malloc(sizeof(float) * _size);
            fftwf_execute_r2r((fftwf_plan)DCT, _weight, cosine);
            double **kernel = new double *[Dimension];
            for (unsigned i = 0; i < Dimension; i++)
            {
                const double pi = 3.14159265358979323846;
                double *k = kernel[i] = new double[points[i]];
                double radius = percent / 100.0f * points[i];
                for (int j = 0; j < points[i]; j++)
                    k[j] = exp(-j * j * radius * radius * pi * pi * 2);
            }
            for (size_t x = 0; x < _size; x++)
            {
                double k = 1.0;
                for (unsigned i = 0; i < Dimension; i++)
                {
                    unsigned j = (x / stride[i]) % points[i];
                    k *= kernel[i][j];
                }
                cosine[x] *= k;
            }
            fftwf_execute_r2r((fftwf_plan)DCT, cosine, _density);
            fftwf_free(cosine);
            for (unsigned i = 0; i < Dimension; i++)
                free(kernel[i]);
            delete[] kernel;
        }
        else
        {
            std::copy(_weight, _weight + _size, _density);
        }

        std::vector<float> sorted;
        std::copy(_density, _density + _size, std::back_inserter(sorted));
        std::sort(sorted.begin(), sorted.end());
        for (size_t x = 0; x < _size; x++)
            _P[x] = (std::lower_bound(sorted.begin(), sorted.end(), _density[x]) - sorted.begin()) / static_cast<float>(_size);

        for_each_fiber([this](Fiber &fiber)
                       { this->basis_functions(fiber); });
        for_each_fiber([this](Fiber &fiber)
                       { this->natural_parameters(fiber); });
        for (size_t x = 0; x < _size; x++)
        {
            for (unsigned j = 0; j < Dimension; j++)
                for (unsigned i = 0; i < Dimension; i++)
                    _T[j][x] += _t[i][j][x];
            for (unsigned j = 0; j < Dimension; j++)
                _L[x] += _T[j][x];
        }
    }

    Covariant(int points[Dimension]) : points(points)
    {
        init();
    }

    Covariant(int grid = 256)
    {
        for (unsigned i = 0; i < Dimension; i++)
            points[i] = grid + 1;
        init();
    }

    ~Covariant()
    {
        fftwf_destroy_plan((fftwf_plan)DCT);
        fftwf_free(_weight);
        fftwf_free(_density);
    };

private:
    void *basis_functions(Fiber &fiber)
    {
        double marginal = (_density[fiber.base] + _density[fiber.base + (points[fiber.d] - 1) * fiber.stride]) / 2.0;
        for (int i = 1; i < points[fiber.d] - 1; i++)
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
            unsigned m = points[fiber.d] / 2;
            for (int j = 0; j < points[fiber.d]; j++)
            {
                if (_f[i][fiber.base + j * fiber.stride] > max)
                {
                    max = _f[i][fiber.base + j * fiber.stride];
                    m = j;
                }
            }
            double delta = 1.0 / (points[fiber.d] - 1);
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
                    t = -2.0 * (std::log(_f[i][fiber.base + j * fiber.stride]) - std::log(_f[i][fiber.base + k * fiber.stride]) / delta / delta / (j - k) / (j - k) - _s[i][fiber.d][fiber.base + k * fiber.stride] / delta / (j - k));
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
                    t = -2.0 * (std::log(_f[i][fiber.base + j * fiber.stride]) - std::log(_f[i][fiber.base + k * fiber.stride]) / delta / delta / (j - k) / (j - k) - _s[i][fiber.d][fiber.base + k * fiber.stride] / delta / (j - k));
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

            for (fiber.id = 0; fiber.id < _size / points[fiber.d]; fiber.id++)
            {
                size_t smaller = dim == 0 ? 1 : fiber.id % stride[dim];
                size_t larger = fiber.id / stride[dim];
                fiber.base = larger * points[fiber.d] * smaller;
                func(fiber);
            }
        }
    }

    void init()
    {
        _size = 1;
        for (unsigned i = 0; i < Dimension; i++)
        {
            stride[i] = _size;
            kind[i] = FFTW_REDFT00;
            _size *= points[i];
        }
        for (unsigned i = 0; i < Dimension; i++)
        {
            _f[i].resize(_size);
            for (unsigned j = 0; j < Dimension; j++)
            {
                _s[i][j].resize(_size);
                _t[i][j].resize(_size);
            }
        }
        _L.resize(_size);
        _P.resize(_size);
        _weight = (float *)fftw_malloc(sizeof(float) * _size);
        _density = (float *)fftw_malloc(sizeof(float) * _size);
        DCT = (void *)fftwf_plan_r2r(Dimension, points, _weight, _density, kind, 0);
        assert(DCT);
    }
};
