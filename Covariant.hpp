#pragma once

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
    int points[Dimension];

private:
    struct Fiber
    {
        Covariant &cov;
        size_t id;
        size_t base;
        size_t stride;
        double delta;
        unsigned d;
        float &f(int i, int j)
        {
            return cov._f[i][base + j * stride];
        }
        float &s(int i, int j)
        {
            return cov._s[i][d][base + j * stride];
        }
        float &t(int i, int j)
        {
            return cov._t[i][d][base + j * stride];
        }
        float &S(int j)
        {
            return cov._S[d][base + j * stride];
        }
        float &T(int j)
        {
            return cov._T[d][base + j * stride];
        }
        float &P(int j)
        {
            return cov._QC[base + j * stride];
        }
        Fiber(Covariant &cov) : cov(cov) {};
    };

    size_t _size;
    size_t _events = 0;
    float *_weight;
    float *_density;
    size_t stride[Dimension];
    fftw_r2r_kind kind[Dimension];
    void *DCT = nullptr;
    unsigned long fft_normalizer = 1;
    std::array<std::vector<float>, Dimension> _f;
    std::array<std::array<std::vector<float>, Dimension>, Dimension> _s;
    std::array<std::array<std::vector<float>, Dimension>, Dimension> _t;
    std::array<std::vector<float>, Dimension> _S;
    std::array<std::vector<float>, Dimension> _T;
    std::vector<Covariant<Dimension>::Event> _M;
    std::vector<float> _L;
    std::vector<float> _QC;
    std::array<std::vector<float>, Dimension> _r;
    std::array<std::vector<float>, Dimension> _q;
    std::array<std::vector<float>, Dimension> _P;
    std::array<std::vector<float>, Dimension> _Q;
    std::vector<float> _R;

    public:
    size_t size() const
    {
        return _size;
    }

    size_t events()
    {
        return _events;
    }

    const float *w() const
    {
        return _weight;
    }

    const float &w(size_t x) const
    {
        return _weight[x];
    }

    const float *f() const
    {
        return _density;
    }

    const float *f(unsigned i) const
    {
        return _f[i].data();
    }

    const float &f(unsigned i, size_t x) const
    {
        assert(x < _size);
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

    const float *S(unsigned i) const
    {
        return _S[i].data();
    }

    const float &S(unsigned i, size_t x) const
    {
        return _S[i][x];
    }

    const float *T(unsigned i) const
    {
        return _T[i].data();
    }

    const float &T(unsigned i, size_t x) const
    {
        return _T[i][x];
    }

    const float &L(size_t x) const
    {
        return _L[x];
    }

    const float * QC() const
    {
        return _QC.data();
    }

    const float &QC(size_t x) const
    {
        return _QC[x];
    }

    const float *R() const
    {
        return _R.data();
    }

    const float &R(size_t x) const
    {
        return _R[x];
    }

    const float *P(unsigned i) const
    {
        return _P[i].data();
    }
    const float *Q(unsigned i) const
    {
        return _Q[i].data();
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
            assert(x + offset < _size);
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
            for (size_t x = 0; x < _size; x++)
                if (_density[x] < 0.0f)
                    _density[x] = 0.0f;
        }
        else
        {
            std::copy(_weight, _weight + _size, _density);
        }

        std::vector<float> sorted;
        std::copy(_density, _density + _size, std::back_inserter(sorted));
        std::sort(sorted.begin(), sorted.end());
        std::vector<float> summed;
        summed.resize(sorted.size());
        double sum = 0.0;
        for (size_t x = 0; x < _size; x++)
            summed[x] = sum += sorted[x];
        for (size_t x = 0; x < _size; x++)
        {
            _QC[x] = summed.at((std::lower_bound(sorted.begin(), sorted.end(), _density[x]) - sorted.begin())) / sum;
            _density[x] /= sum;
            size_t y = x;
            for (unsigned i = 0, j; i < Dimension; i++)
            {
                j = y % points[i];
                y /= points[i];
                _P[i][j] += _density[x];
            }
        }

        for_each_fiber([this](Fiber &fiber)
                       { this->basis_functions(fiber); });
        for_each_fiber([this](Fiber &fiber)
                       { this->natural_parameters(fiber); });

        for (size_t x = 0; x < _size; x++)
        {
            size_t y = x;
            for (unsigned i = 0, j; i < Dimension; i++)
            {
                double dual = 1.0;
                for (j = 0; j < Dimension; j++)
                {
                    if (j == i)
                        continue;
                    if (_f[i][x] > 0.0f)
                        _r[i][x] += (_t[i][j][x] + _t[j][i][x]) / _f[i][x] / _f[i][x];
                    for (unsigned k = 0; k < Dimension; k++)
                    {
                        if (k == i || k == j)
                            continue;
                        if (_f[j][x] > 0.0f)
                            _q[i][x] += (_t[j][k][x] + _t[k][j][x]) / _f[j][x] / _f[j][x];
                    }
                    dual *= _f[j][x];
                }
                j = y % points[i];
                y /= points[i];
                _R[x] += _r[i][x] * _density[x];
                _Q[i][j] += _q[i][x] * dual * _P[i][j];
            }
        }

        for (size_t x = 0; x < _size; x++)
            for (unsigned i = 0; i < Dimension; i++)
                for (unsigned j = 0; j < Dimension; j++)
                {
                    _S[j][x] += _s[i][j][x];
                    _T[j][x] += _s[i][j][x];
                }
        for (size_t x = 0; x < _size; x++)
            for (unsigned j = 0; j < Dimension; j++)
            {
                _L[x] += _T[j][x];
                if (_QC[x] <= 0.05f)
                    _L[x] = -1.0f;
            }
        for (unsigned j = 0; j < Dimension; j++)
        {
            filter(_S[j].data(), percent, true);
            filter(_T[j].data(), percent, true);
        }
    }

    double factorProbability()
    {
        double error = 0.0;
        for (size_t x = 0; x < _size; x++)
        {
            double product = 1.0;
            for (unsigned i = 0; i < Dimension; i++)
                product *= f(i, x);
            double diff = product - _density[x];
            diff = std::abs(diff);
            error += diff * _density[x];
        }
        return error;
    }

    double differentialEquation()
    {
        double error = 0.0;
        for_each_fiber([this, &error](Fiber &fiber)
                       {
            double t, s_diff, t_diff, delta = 1.0 / (points[fiber.d] - 1);
            for (unsigned int i = 0; i < Dimension; i++)
            {
                for (int k = 0; k < points[fiber.d] - 1; k++)
                {
                    if (fiber.f(i, k + 1) > 0 && fiber.f(i, k) > 0) {
                        t = -2.0 * ((std::log(fiber.f(i, k + 1)) - std::log(fiber.f(i, k))) / delta / delta - fiber.s(i, k) / delta);
                        t_diff = std::abs((t - fiber.t(i, k))/fiber.t(i,k));
                    }
                    s_diff = std::abs((fiber.s(i, k + 1) - (-fiber.t(i, k)/(points[fiber.d] - 1) + fiber.s(i, k))) / fiber.s(i, k));
                    if (s_diff > error)
                        error = s_diff;
                    if (t_diff > error)
                        error = t_diff;
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
    void filter(float *input, float *output, float percent = 1.0f, bool normalize = false)
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
        if (normalize)
            for (unsigned x = 0; x < _size; x++)
                output[x] /= (double)fft_normalizer;
    }

    void filter(float *data, float percent = 1.0f, bool normalize = false)
    {
        filter(data, data, percent, normalize);
    }

    void *basis_functions(Fiber &fiber)
    {
        double marginal = (_density[fiber.base] + _density[fiber.base + (points[fiber.d] - 1) * fiber.stride]) / 2.0;
        for (int i = 1; i < points[fiber.d] - 1; i++)
        {
            marginal += _density[fiber.base + i * fiber.stride];
        }
        if (marginal < 1.0 / (double)_events)
            for (int i = 0; i < points[fiber.d]; i++)
                fiber.f(fiber.d, i) = 0.0f;
        else
            for (int i = 0; i < points[fiber.d]; i++)
                fiber.f(fiber.d, i) = _density[fiber.base + i * fiber.stride] / marginal;
        return nullptr;
    }

    inline double squared(double x) { return x * x; };

    void *natural_parameters(Fiber &fiber)
    {
        for (unsigned i = 0; i < Dimension; i++)
        {
            float max = 0.0f;
            if (fiber.d > 0)
                max = 0.0f;
            int m = points[fiber.d] - 1;
            for (int k = 0; k < points[fiber.d]; k++)
            {
                if (fiber.f(i, k) > max)
                {
                    max = fiber.f(i, k);
                    m = k;
                }
            }
            double t;
            fiber.s(i, m) = 0.0;
            for (int k = m, j; k < points[fiber.d] - 1;)
            {
                j = k + 1;
                while (fiber.f(i, j) <= 0.0 && j < points[fiber.d] - 1)
                    j++;
                if (fiber.f(i, j) <= 0.0)
                    t = 1.0 / squared(fiber.delta);
                else
                    t = -2.0 * ((std::log(fiber.f(i, j)) - std::log(fiber.f(i, k))) / squared(fiber.delta * (j - k)) - fiber.s(i, k) / fiber.delta / (j - k));
                while (k < j)
                {
                    fiber.t(i, k) = (float)t;
                    if (k != points[fiber.d] - 1)
                    {
                        double s = -t * fiber.delta + fiber.s(i, k);
                        fiber.s(i, k + 1) = (float)s;
                    }
                    k++;
                }
            }
            for (unsigned k = m, j; k > 0;)
            {
                j = k - 1;
                while (j > 0 && fiber.f(i, j) <= 0.0)
                    j--;
                if (fiber.f(i, j) <= 0.0)
                    t = 1.0 / squared(fiber.delta);
                else
                    t = 2.0 * ((std::log(fiber.f(i, k)) - std::log(fiber.f(i, j))) / squared(fiber.delta * (k - j)) - fiber.s(i, k) / fiber.delta / (k - j));
                while (k > j)
                {
                    fiber.t(i, k - 1) = (float)t;
                    double s = t * fiber.delta + fiber.s(i, k);
                    fiber.s(i, k - 1) = (float)s;
                    k--;
                }
            }
        }
        return nullptr;
    }

    void for_each_fiber(std::function<void(Fiber &)> func)
    {
        Fiber fiber(*this);
        for (unsigned i = 0; i < Dimension; ++i)
        {
            fiber.d = i;
            fiber.stride = stride[i];
            fiber.delta = 1.0 / (double)(points[fiber.d] - 1);

            for (fiber.id = 0; fiber.id < _size / points[fiber.d]; fiber.id++)
            {
                size_t smaller = fiber.id % stride[fiber.d];
                size_t larger = fiber.id / stride[fiber.d];
                if (fiber.d == Dimension - 1)
                    fiber.base = smaller;
                else
                    fiber.base = larger * points[fiber.d] * stride[fiber.d] + smaller;
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
            fft_normalizer *= 2 * (points[i] - 1);
            _f[i].resize(_size);
            _T[i].resize(_size);
            _S[i].resize(_size);
            _r[i].resize(_size);
            _q[i].resize(_size);
            _P[i].resize( points[i] );
            _Q[i].resize( points[i] );
            for (unsigned j = 0; j < Dimension; j++)
            {
                _s[i][j].resize(_size);
                _t[i][j].resize(_size);
            }
        }
        _L.resize(_size);
        _QC.resize(_size);
        _R.resize(_size);
        _weight = (float *)fftw_malloc(sizeof(float) * _size);
        _density = (float *)fftw_malloc(sizeof(float) * _size);
        DCT = (void *)fftwf_plan_r2r(Dimension, (const int *)&points, _weight, _density, kind, 0);
        assert(DCT);
    }
};
