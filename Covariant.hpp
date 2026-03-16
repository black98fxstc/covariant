#include <cmath>
#include <array>
#include <vector>
#include <fftw3.h>
#include <assert.h>

template <unsigned Dimension>
class Covariant
{
public:
    typedef std::array<float, Dimension> Event;
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
    std::array<std::vector<float>, Dimension> _S;
    std::array<std::vector<float>, Dimension> _T;
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

    const float &w(size_t x) const
    {
        return _weight[x];
    }

    const float &f(size_t x) const
    {
        return _density[x];
    }

    const float &f(unsigned i, size_t x) const
    {
        return _f[i][x];
    }

    const float &s(unsigned i, unsigned j, size_t x) const
    {
        return _s[i][j][x];
    }

    const float &t(unsigned i, unsigned j, size_t x) const
    {
        return _t[i][j][x];
    }

    const float &S(unsigned i, size_t x) const
    {
        return _S[i][x];
    }

    const float &T(unsigned i, size_t x) const
    {
        return _T[i][x];
    }

    const float &L(size_t x) const
    {
        return _L[x];
    }

    const float &P(size_t x) const
    {
        return _P[x];
    }

    bool event(const Event &event)
    {
        size_t x = 0;
        double rem[Dimension];
        for (unsigned i = 0; i < Dimension; i++)
        {
            if (event[i] < 0.0f || event[i] >= 1.0f)
                return false;
            rem[i] = event[i] * (points[i] - 1);
            unsigned floor = static_cast<unsigned>(rem[i]);
            rem[i] -= floor;
            x += floor * stride[i];
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
            filter(_weight, _density, percent = 1.0f);
        }
        else
        {
            std::copy(_weight, _weight + _size, _density);
        }

        std::vector<float> sorted;
        std::copy(_density, _density + _size, std::back_inserter(sorted));
        std::sort(sorted.begin(), sorted.end());
        std::vector<float> summed(sorted);
        double sum = 0.0;
        for (size_t x = 0; x < _size; x++)
            summed[x] = sum += sorted[x];
        for (size_t x = 0; x < _size; x++)
        {
            _P[x] = summed[(std::lower_bound(sorted.begin(), sorted.end(), _density[x]) - sorted.begin())] / sum;
            if (_density[x] < 0.0f)
                _density[x] = 0.0f;
            else
                _density[x] /= sum;
        }

        for_each_fiber([this](Fiber &fiber)
                       { this->basis_functions(fiber); });
        for_each_fiber([this](Fiber &fiber)
                       { this->natural_parameters(fiber); });
        for (size_t x = 0; x < _size; x++)
        {
            for (unsigned j = 0; j < Dimension; j++) {
                for (unsigned i = 0; i < Dimension; i++)
                {
                    _S[j][x] += _s[i][j][x];
                    _T[j][x] += _t[i][j][x];
                }
                if (_P[x] <= 0.05f)
                    _S[j][x] = 1.0f;
            }
            for (unsigned j = 0; j < Dimension; j++)
                _L[x] += _T[j][x];
        }

        // for_each_fiber([this](Fiber &fiber)
        // {
        //     const double threshold = 0.05f;
        //     int p = 0, q;
        //     for (int i = 0; i < points[fiber.d]; i++)
        //         if (_P[fiber.base + i * fiber.stride] <= threshold)
        //         _S[fiber.d][fiber.base + i * fiber.stride] = -1.0f;
        //     // while (_P[fiber.base + p * fiber.stride] <= threshold && p < points[fiber.d])
        //     //     p++;
        //     // if (p == points[fiber.d]) {
        //     //     for (int i = 0; i < points[fiber.d]; i++)
        //     //         _S[fiber.d][fiber.base + i * fiber.stride] = -1.0f;
        //     //     return;     
        //     // }
        //     // float edge_value = _S[fiber.d][fiber.base + p * fiber.stride];
        //     // for (int i = p; i > 0; i--)
        //     //     _S[fiber.d][fiber.base + i * fiber.stride] = edge_value;
        //     // for (q = p + 1; p < points[fiber.d]; p = q) {
        //     //     while (_P[fiber.base + q * fiber.stride] > threshold && q < points[fiber.d])
        //     //         q++;
        //     //     if (q == points[fiber.d])
        //     //         break;
        //     //     p = q;
        //     //     while (_P[fiber.base + q * fiber.stride] <= threshold && q < points[fiber.d])
        //     //         q++;
        //     //     if (q == points[fiber.d])
        //     //         break;
        //     //     for (int i = p + 1; i < q; i++)
        //     //         _S[fiber.d][fiber.base + i * fiber.stride] =
        //     //             (double)(i-p)/(double)(q-p)*_S[fiber.d][fiber.base + p * fiber.stride]
        //     //           + (double)(q-i)/(double)(q-p)*_S[fiber.d][fiber.base + q * fiber.stride];
        //     // }
        //     // edge_value = _S[fiber.d][fiber.base + p * fiber.stride];
        //     // for (int i = p; i < points[fiber.d]; i++)
        //     //     _S[fiber.d][fiber.base + i * fiber.stride] = edge_value;
        //     });
    }

    double factorProbability()
    {
        double error = 0.0;
        for (size_t x = 0; x < _size; x++)
        {
            double product = 1.0;
            for (unsigned i = 0; i < Dimension; i++)
                product *= f(i, x);
            double diff = product - f(x);
            diff = std::abs(diff);
            if (diff > error)
                error = diff;
        }
        return error;
    }

    double differentialEquation()
    {
        double error = 0.0;
        for_each_fiber([this, &error](Fiber &fiber)
                       {
            for (unsigned i = 0; i < Dimension; i++)
            {
                for (int k = 0; k < points[fiber.d] - 1; k++)
                {
                    double diff = _s[i][fiber.d][fiber.base + (k + 1) * fiber.stride] - (-_t[i][fiber.d][fiber.base + k * fiber.stride] * 1.0/(points[fiber.d] - 1) + _s[i][fiber.d][fiber.base + k * fiber.stride]);
                    diff = std::abs(diff);
                    if (diff > error)
                        error = diff;
                }
            } });
        return error;
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
    }

private:
    void filter(float* input, float* output, float percent = 1.0f)
     {
        float *cosine = (float *)fftw_malloc(sizeof(float) * _size);
        fftwf_execute_r2r((fftwf_plan)DCT, input, cosine);
        double **kernel = new double *[Dimension];
        for (unsigned i = 0; i < Dimension; i++)
        {
            const double pi = 3.14159265358979323846;
            double *k = kernel[i] = new double[points[i]];
            double radius = percent / 100.0f;
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
        fftwf_execute_r2r((fftwf_plan)DCT, cosine, output);
        fftwf_free(cosine);
        for (unsigned i = 0; i < Dimension; i++)
            free(kernel[i]);
        delete[] kernel;
    }

    void filter(float* data, float percent = 1.0f)
    {
        filter(data, data, percent);
    }

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
            unsigned m = points[fiber.d];
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
                if (_f[i][fiber.base + j * fiber.stride] <= 0.0)
                    t = 1.0; // delta / delta;
                else
                    t = -2.0 * ((std::log(_f[i][fiber.base + j * fiber.stride]) - std::log(_f[i][fiber.base + k * fiber.stride])) / delta / delta / (j - k) / (j - k) - _s[i][fiber.d][fiber.base + k * fiber.stride] / delta / (j - k));
                while (k < j)
                {
                    _t[i][fiber.d][fiber.base + k * fiber.stride] = t;
                    if (k != points[fiber.d] - 1)
                        _s[i][fiber.d][fiber.base + (k + 1) * fiber.stride] = -t * delta + _s[i][fiber.d][fiber.base + k * fiber.stride];
                    k++;
                }
            }
            for (unsigned k = m, j; k > 0;)
            {
                j = k - 1;
                while (_f[i][fiber.base + j * fiber.stride] <= 0.0 && j > 0)
                    j--;
                if (_f[i][fiber.base + j * fiber.stride] <= 0.0)
                    t = 1.0 / delta / delta;
                else
                    t = -2.0 * ((std::log(_f[i][fiber.base + j * fiber.stride]) - std::log(_f[i][fiber.base + k * fiber.stride])) / delta / delta / (j - k) / (j - k) - _s[i][fiber.d][fiber.base + k * fiber.stride] / delta / (j - k));
                while (k > j)
                {
                    _t[i][fiber.d][fiber.base + (k - 1) * fiber.stride] = t;
                    _s[i][fiber.d][fiber.base + (k - 1) * fiber.stride] = t * delta + _s[i][fiber.d][fiber.base + k * fiber.stride];
                    k--;
                }
            }
        }
        return nullptr;
    }

    void for_each_fiber(std::function<void(Fiber &)> func)
    {
        Fiber fiber;
        for (unsigned i = 0; i < Dimension; ++i)
        {
            fiber.d = i;
            fiber.stride = stride[i];

            for (fiber.id = 0; fiber.id < _size / points[fiber.d]; fiber.id++)
            {
                size_t smaller = fiber.id % stride[i];
                size_t larger = fiber.id / stride[i];
                fiber.base = larger * stride[i] * points[fiber.d] + smaller;
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
            _T[i].resize(_size);
            _S[i].resize(_size);
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
