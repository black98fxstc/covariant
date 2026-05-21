#pragma once

#include <cmath>
#include <array>
#include <vector>
#include <fftw3.h>
#include <assert.h>
#include <type_traits>

template <unsigned Dimension, typename Float = float>
class Covariant
{
public:
    typedef std::array<Float, Dimension> Event;
    const unsigned dimension = Dimension;
    unsigned points[Dimension];

private:
    struct Fiber
    {
        Covariant<Dimension, Float> &cov;
        size_t id;
        size_t base;
        size_t stride;
        double delta;
        unsigned d;
        Float &f(int i, int j)
        {
            return cov._f[i][base + j * stride];
        }
        Float &s(int i, int j)
        {
            return cov._s[i][d][base + j * stride];
        }
        Float &t(int i, int j)
        {
            return cov._t[i][d][base + j * stride];
        }
        Float &S(int j)
        {
            return cov._S[d][base + j * stride];
        }
        Float &T(int j)
        {
            return cov._T[d][base + j * stride];
        }
        Float &P(int j)
        {
            return cov._QC[base + j * stride];
        }
        Fiber(Covariant &cov) : cov(cov) {};
    };

    class Hypercube
    {
    private:
        Covariant<Dimension, Float> &cov;

    public:
        enum state { unknown, outlier, ambiguous, contiguous, assigned };

        size_t x;

        Hypercube(Covariant<Dimension, Float> &cov, size_t x) : cov(cov), x(x) 
        {
            cov._status[x] = unknown;
        };

        Float density() const
        {
            return cov._density[x];
        }

        Float quantile() const
        {
            return cov._QC[x];
        }

        Float lapacian() const
        {
            return cov._L[x];
        }

        Float &cluster()
        {
            return cov._cluster[x];
        }

        state &status() const
        {
            return cov._status[x];
        }

        Hypercube &operator=(const Hypercube &other)
        {
            x = other.x;
            return *this;
        }
    };

    size_t _size;
    size_t _events = 0;
    Float *_weight;
    Float *_density;
    size_t stride[Dimension];
    fftw_r2r_kind kind[Dimension];
    void *DCT = nullptr;
    Float _tot_R = 0.0f;
    Float s_max = std::numeric_limits<Float>::lowest(), s_min = std::numeric_limits<Float>::max();
    Float t_max = std::numeric_limits<Float>::lowest(), t_min = std::numeric_limits<Float>::max();
    double _differential_error = 0.0;
    double _factor_error = 0.0;
    unsigned long fft_normalizer = 1;
    std::vector<Float> _L;
    std::vector<Float> _QC;
    std::vector<Float> _R;
    std::vector<Float> _cluster;
    std::vector<typename Hypercube::state> _status;
    std::array<std::vector<Float>, Dimension> _f;
    std::array<std::vector<Float>, Dimension> _S;
    std::array<std::vector<Float>, Dimension> _T;
    std::array<std::vector<Float>, Dimension> _r;
    std::array<std::vector<Float>, Dimension> _q;
    std::array<std::vector<Float>, Dimension> _P;
    std::array<std::vector<Float>, Dimension> _Q;
    std::array<std::vector<Float>, Dimension> _var_Q;
    std::array<Float, Dimension> _tot_Q;
    std::array<Float, Dimension> _var_tot_Q;
    std::array<std::array<Float*, Dimension>, Dimension> _s;
    std::array<std::array<Float*, Dimension>, Dimension> _t;

public:
    size_t size() const
    {
        return _size;
    }

    size_t events()
    {
        return _events;
    }

    const Float *w() const
    {
        return _weight;
    }

    const Float &w(size_t x) const
    {
        return _weight[x];
    }

    const Float *f() const
    {
        return _density;
    }

    const Float *f(unsigned i) const
    {
        return _f[i].data();
    }

    const Float &f(unsigned i, size_t x) const
    {
        assert(x < _size);
        return _f[i][x];
    }

    const Float &s(unsigned i, unsigned j, size_t x) const
    {
        return _s[i][j][x];
    }

    const Float *t(unsigned i, unsigned j) const
    {
        return _t[i][j];
    }

    const Float &t(unsigned i, unsigned j, size_t x) const
    {
        return _t[i][j][x];
    }

    const Float *S(unsigned i) const
    {
        return _S[i].data();
    }

    const Float &S(unsigned i, size_t x) const
    {
        return _S[i][x];
    }

    const Float *T(unsigned i) const
    {
        return _T[i].data();
    }

    const Float &T(unsigned i, size_t x) const
    {
        return _T[i][x];
    }

    const Float *L() const
    {
        return _L.data();
    }

    const Float &L(size_t x) const
    {
        return _L[x];
    }

    const Float *QC() const
    {
        return _QC.data();
    }

    const Float &QC(size_t x) const
    {
        return _QC[x];
    }

    const Float *R() const
    {
        return _R.data();
    }

    const Float &R(size_t x) const
    {
        return _R[x];
    }

    const Float &tot_R() const
    {
        return _tot_R;
    }

    const Float *P(unsigned i) const
    {
        return _P[i].data();
    }
    const Float *Q(unsigned i) const
    {
        return _Q[i].data();
    }

    const Float *classes() const
    {
        return _cluster.data();
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

    unsigned short classify(const Event &event)
    {
        size_t x = 0;
        for (unsigned i = 0; i < Dimension; i++)
        {
            if (event[i] < 0.0f || event[i] >= 1.0f)
                return 0;
            unsigned floor = static_cast<unsigned>(event[i] * (points[i] - 1));
            x += floor * stride[i];
        }
        return (int)(_cluster[x]);
    }

    void parameters(Float smoothing = .01f, Float threshold = 0.001f)
    {
        if (smoothing > 0.0f)
        {
            filter(_weight, _density, smoothing);
            for (size_t x = 0; x < _size; x++)
                if (_density[x] < 0.0f)
                    _density[x] = 0.0f;
        }
        else
        {
            std::copy(_weight, _weight + _size, _density);
        }

        std::vector<Float> sorted;
        std::copy(_density, _density + _size, std::back_inserter(sorted));
        std::sort(sorted.begin(), sorted.end());
        std::vector<Float> summed;
        summed.resize(sorted.size());
        double sum = 0.0;
        for (size_t x = 0; x < _size; x++)
            summed[x] = sum += sorted[x];
        for (size_t x = 0; x < _size; x++)
        {
            _QC[x] = summed.at((std::lower_bound(sorted.begin(), sorted.end(), _density[x]) - sorted.begin())) / sum;
            _density[x] /= sum;
            for (unsigned i = 0; i < Dimension; i++)
            {
                unsigned j = (x / stride[i]) % points[i];
                _P[i][j] += _density[x];
            }
        }

        for_each_fiber([this](Fiber &fiber)
                       { this->basis_functions(fiber); });
        for_each_fiber([this](Fiber &fiber)
                       { this->natural_parameters(fiber); });

        for (size_t x = 0; x < _size; x++)
        {
            if (_density[x] <= 0.0f)
                continue;
            double product = 1.0;
            for (unsigned i = 0; i < Dimension; i++)
                product *= _f[i][x];
            double error = std::abs(product - _density[x]);
            _factor_error += error * _density[x];
        }
        for_each_fiber([this](Fiber &fiber)
                       { this->differential_error(fiber); });

        if (smoothing > 0.0f)
            for (unsigned j = 0; j < Dimension; j++)
                for (unsigned i = 0; i < Dimension; i++)
                {
                    filter(_s[i][j], 1.0f/(float)(points[j] - 1), true);
                    filter(_t[i][j], 1.0f/(float)(points[j] - 1), true);
                }

        for (size_t x = 0; x < _size; x++)
        {
            double dual;
            for (unsigned i = 0; i < Dimension; i++)
            {
                dual = 1.0;
                for (unsigned j = 0; j < Dimension; j++)
                {
                    _S[j][x] += _s[i][j][x];
                    _T[j][x] += _t[i][j][x];

                    if (j == i)
                        continue;
                    if (_f[i][x] > 0.0f)
                        _r[i][x] += (_t[i][j][x] + _t[j][i][x]);// / squared(_f[i][x]);
                    for (unsigned k = 0; k < Dimension; k++)
                    {
                        if (k == i || k == j)
                            continue;
                        if (_f[j][x] > 0.0f)
                            _q[i][x] += (_t[j][k][x] + _t[k][j][x]);// / squared(_f[j][x]);
                    }
                    dual *= _f[j][x];
                }
                _R[x] += _r[i][x];

                if (_QC[x] < threshold)
                    continue;
                unsigned j = (x / stride[i]) % points[i];
                _Q[i][j] += _q[i][x] * dual;
                _tot_Q[i] += _Q[i][j] * _P[i][j];
            }
            for (unsigned i = 0; i < Dimension; i++)
                _L[x] += _T[i][x];
            _R[x] *= 0.25f;// _density[x] / (Float)_size;
            if (_QC[x] >= threshold)
                _tot_R += _R[x];
        }

        for (size_t x = 0; x < _size; x++)
        {
            if (_QC[x] < threshold)
                continue;
            for (unsigned i = 0; i < Dimension; i++)
            {
                unsigned j = (x / stride[i]) % points[i];
                _var_Q[i][j] += squared(_Q[i][j] - _tot_Q[i]);
                _var_tot_Q[i] = _var_Q[i][j] * _P[i][j];
            }
        }

        trim(_L, threshold);
        trim(_R, threshold);
        for (unsigned i = 0; i < Dimension; i++)
        {
            trim(_f[i], threshold);
            trim(_S[i], threshold);
            trim(_T[i], threshold);
            for (unsigned j = 0; j < Dimension; j++)
            {
                trim(_s[i][j], threshold);
                trim(_t[i][j], threshold);
            }
        }
    }

    unsigned bluster(Float threshold = 0.001f)
    {
        for (unsigned x = 0; x < _size; x++)
            if (std::isnan(_L[x]))
                _cluster[x] = std::numeric_limits<Float>::quiet_NaN();
            else if (_L[x] <= 0.0f)
                _cluster[x] = std::numeric_limits<Float>::quiet_NaN();
            else
                _cluster[x] = _L[x];
        return 0;
    }

    unsigned cluster(Float threshold = 0.001f, bool grow = false)
    {
        unsigned clusters = 0;
        std::vector<Hypercube> cubes;
        cubes.reserve(_size);
        for (size_t x = 0; x < _size; x++)
            cubes.emplace_back(*this, x);
        auto outliers = std::partition(cubes.begin(), cubes.end(), [threshold](const Hypercube &cube)
                       { return cube.quantile() >= threshold; });
        auto precision = std::partition(cubes.begin(), outliers, [](const Hypercube &cube)
                       { return cube.lapacian() > 0.0f; });
        std::sort(cubes.begin(), precision, [](const Hypercube &a, const Hypercube &b)
                  { return a.density() > b.density(); });
        for (auto clustered = cubes.begin(); clustered != precision;)
        {
            clusters++;
            clustered->status() = Hypercube::contiguous;
            auto contiguous = clustered + 1;
            while (true)
            {
                for (auto it = clustered; it != contiguous; ++it)
                {
                    it->status() = Hypercube::assigned;
                    it->cluster() = clusters;
                    size_t x = it->x, y = it->x;
                    for (unsigned i = 0; i < Dimension; i++)
                    { 
                        unsigned j = y % points[i];
                        y /= points[i];
                        if (j < points[i] - 1 && _status[x + stride[i]] == Hypercube::unknown)
                            _status[x + stride[i]] = Hypercube::contiguous;
                        if (j > 0 && _status[x - stride[i]] == Hypercube::unknown)
                            _status[x - stride[i]] = Hypercube::contiguous;
                    }
                }
                clustered = contiguous;
                contiguous = std::partition(clustered, precision, [](const Hypercube &cube)
                        { return cube.status() == Hypercube::contiguous; });
                if (contiguous == clustered)
                    break;
            }
        }
        if (!grow)
        {
            for (auto it = precision; it != outliers; ++it)
            {
                it->status() = Hypercube::ambiguous;
                it->cluster() = std::numeric_limits<Float>::quiet_NaN();
            }
            for (auto it = outliers; it != cubes.end(); ++it)
                {
                    it->status() = Hypercube::outlier;
                    it->cluster() = std::numeric_limits<Float>::quiet_NaN();
                }
        }
        else
        {
            auto contiguous = precision;
            for (auto clustered = cubes.begin(); clustered != cubes.end(); )
            {
                for (auto it = clustered; it != contiguous; ++it)
                {
                    size_t x = it->x, y = it->x;
                    for (unsigned i = 0; i < Dimension; i++)
                    { 
                        unsigned j = y % points[i];
                        y /= points[i];
                        if (j < points[i] - 1 && _status[x + stride[i]] == Hypercube::unknown)
                        {
                            _status[x + stride[i]] = Hypercube::contiguous;
                            _cluster[x + stride[i]] = it->cluster();
                        }
                        if (j > 0 && _status[x + stride[i]] == Hypercube::unknown)
                        {
                            _status[x - stride[i]] = Hypercube::contiguous;
                            _cluster[x - stride[i]] = it->cluster();
                        }
                    }
                }
                clustered = contiguous;
                contiguous = std::partition(clustered, cubes.end(), [](const Hypercube &cube)
                        { return cube.status() == Hypercube::contiguous; });
                if (contiguous == clustered)
                    break;
            }
        }
        return clusters;
    }

    double factorProbability()
    {
        return _factor_error;
    }

    double differentialEquation()
    {
        return _differential_error;
    }

    Covariant(unsigned points[Dimension], bool column_major = false) : points(points)
    {
        init(column_major);
    }

    Covariant(unsigned grid = 256, bool column_major = false)
    {
        for (unsigned i = 0; i < Dimension; i++)
            points[i] = grid + 1;
        init(column_major);
    }

    ~Covariant()
    {
        if constexpr (std::is_same_v<Float, double>)
            fftw_destroy_plan((fftw_plan)DCT);
        else
            fftwf_destroy_plan((fftwf_plan)DCT);
        fftw_free(_weight);
        fftw_free(_density);
        for (unsigned i = 0; i < Dimension; i++)
            for (unsigned j = 0; j < Dimension; j++)
            {
                fftw_free(_s[i][j]);
                fftw_free(_t[i][j]);
            }
    }

private:
    void basis_functions(Fiber &fiber)
    {
        double marginal = (_density[fiber.base] + _density[fiber.base + (points[fiber.d] - 1) * fiber.stride]) / 2.0;
        for (unsigned i = 1; i < points[fiber.d] - 1; i++)
        {
            marginal += _density[fiber.base + i * fiber.stride];
        }
        if (marginal < 1.0 / (double)_events)
            for (unsigned i = 0; i < points[fiber.d]; i++)
                fiber.f(fiber.d, i) = 0.0f;
        else
            for (unsigned i = 0; i < points[fiber.d]; i++)
                fiber.f(fiber.d, i) = _density[fiber.base + i * fiber.stride] / marginal;
    }

    void natural_parameters(Fiber &fiber)
    {
        for (unsigned i = 0; i < Dimension; i++)
        {
            Float max;
            max = 0.0f;
            int m = points[fiber.d] - 1;
            fiber.t(i, m) = 1.0 / fiber.delta;
            for (unsigned k = 0; k < points[fiber.d]; k++)
            {
                if (fiber.f(i, k) > max)
                {
                    max = fiber.f(i, k);
                    m = k;
                }
            }

            double t;
            fiber.s(i, m) = 0.0;
            for (unsigned k = m, j; k < points[fiber.d] - 1;)
            {
                j = k + 1;
                while (fiber.f(i, j) <= 0.0 && j < points[fiber.d] - 1)
                    j++;
                if (fiber.f(i, j) <= 0.0)
                    t = 1.0 / fiber.delta;
                else
                    t = -2.0 * ((std::log(fiber.f(i, j)) - std::log(fiber.f(i, k))) / squared(fiber.delta * (j - k)) - fiber.s(i, k) / fiber.delta / (j - k));
                while (k < j)
                {
                    fiber.t(i, k) = (Float)t;
                    if (t > t_max)
                        t_max = t;
                    if (t < t_min)
                        t_min = t;
                    if (k != points[fiber.d] - 1)
                    {
                        double s = -t * fiber.delta + fiber.s(i, k);
                        if (s > s_max)
                            s_max = s;
                        if (s < s_min)
                            s_min = s;
                        fiber.s(i, k + 1) = (Float)s;
                    }
                    k++;
                }
            }
            for (unsigned k = m, j; k > 0;)
            {
                j = k - 1;
                while (j > 0 && fiber.f(i, j) <= 0.0)
                    j--;
                if (fiber.f(i, j) <= 0.0 || fiber.f(i, k) <= 0.0)
                    t = 1.0 / squared(fiber.delta);
                else
                    t = 2.0 * ((std::log(fiber.f(i, k)) - std::log(fiber.f(i, j))) / squared(fiber.delta * (k - j)) - fiber.s(i, k) / fiber.delta / (k - j));
                if (t > t_max)
                    t_max = t;
                if (t < t_min)
                    t_min = t;
                while (k > j)
                {
                    fiber.t(i, k - 1) = (Float)t;
                    double s = t * fiber.delta + fiber.s(i, k);
                    if (s > s_max)
                        s_max = s;
                    if (s < s_min)
                        s_min = s;
                    fiber.s(i, k - 1) = (Float)s;
                    k--;
                }
            }
        }
    }

    void differential_error(Fiber &fiber)
    {
        double t, s_diff, t_diff;
        for (unsigned int i = 0; i < Dimension; i++)
        {
            for (unsigned k = 0; k < points[fiber.d] - 1; k++)
            {
                if (fiber.f(i, k + 1) > 0 && fiber.f(i, k) > 0) 
                {
                    t = -2.0 * ((std::log(fiber.f(i, k + 1)) - std::log(fiber.f(i, k))) / squared(fiber.delta) - fiber.s(i, k) / fiber.delta);
                    t_diff = std::abs((t - fiber.t(i, k)) / (t_max - t_min));
                    if (t_diff > _differential_error)
                        _differential_error = t_diff;
                }
                if (fiber.s(i, k + 1) != 0)
                {
                    s_diff = std::abs((fiber.s(i, k + 1) - (-fiber.t(i, k) * fiber.delta + fiber.s(i, k))) / (s_max - s_min));
                    if (s_diff > _differential_error)
                        _differential_error = s_diff;
                }
            }
        }
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
                fiber.base = larger * points[fiber.d] * stride[fiber.d] + smaller;
                func(fiber);
            }
        }
    }

    inline double squared(double x) { return x * x; };

    void filter(Float *input, Float *output, Float radius = 0.01f, bool normalize = false)
    {
        Float *cosine = (Float *)fftw_malloc(sizeof(Float) * _size);
        if constexpr (std::is_same_v<Float, double>)
            fftw_execute_r2r((fftw_plan)DCT, input, cosine);
        else
            fftwf_execute_r2r((fftwf_plan)DCT, input, cosine);
        double **kernel = new double *[Dimension];
        for (unsigned i = 0; i < Dimension; i++)
        {
            const double pi = 3.14159265358979323846;
            double *k = kernel[i] = new double[points[i]];
            for (unsigned j = 0; j < points[i]; j++)
                k[j] = exp(-2.0 * squared(j * radius * pi));
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
        if constexpr (std::is_same_v<Float, double>)
            fftw_execute_r2r((fftw_plan)DCT, cosine, output);
        else
            fftwf_execute_r2r((fftwf_plan)DCT, cosine, output);
        fftw_free(cosine);
        for (unsigned i = 0; i < Dimension; i++)
            free(kernel[i]);
        delete[] kernel;
        if (normalize)
            for (unsigned x = 0; x < _size; x++)
                output[x] /= (Float)fft_normalizer;
    }

    void filter(Float *data, Float radius = 0.01f, bool normalize = false)
    {
        filter(data, data, radius, normalize);
    }

    void trim(std::vector<Float> &data, Float threshold)
    {
        trim(data.data(), threshold);
    }

    void trim(Float *data, Float threshold)
    {
        for (size_t x = 0; x < _size; x++)
            if (_QC[x] < threshold)
                data[x] = std::numeric_limits<Float>::quiet_NaN();
    }

    void init(bool column_major)
    {
        _size = 1;
        if (column_major)
            for (int i = Dimension - 1; i >= 0; i--)
            {
                stride[i] = _size;
                kind[i] = FFTW_REDFT00;
                _size *= points[i];
            }
        else
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
            _P[i].resize(points[i]);
            _Q[i].resize(points[i]);
            _var_Q[i].resize(points[i]);
            for (unsigned j = 0; j < Dimension; j++)
            {
                _s[i][j] = (Float *)fftw_malloc(sizeof(Float) * _size);
                _t[i][j] = (Float *)fftw_malloc(sizeof(Float) * _size);
            }
        }
        _L.resize(_size);
        _QC.resize(_size);
        _R.resize(_size);
        _cluster.resize(_size);
        _status.resize(_size);
        _weight = (Float *)fftw_malloc(sizeof(Float) * _size);
        _density = (Float *)fftw_malloc(sizeof(Float) * _size);
        assert((std::is_same_v<Float, float> || std::is_same_v<Float, double>));
        int fftw_n[Dimension];
        std::copy(points, points + Dimension, fftw_n);
        if (column_major)
            std::reverse(fftw_n, fftw_n + Dimension);
        if constexpr (std::is_same_v<Float, double>)
            DCT = (void *)fftw_plan_r2r(Dimension, fftw_n, _weight, _density, kind, 0);
        else
            DCT = (void *)fftwf_plan_r2r(Dimension, fftw_n, _weight, _density, kind, 0);
        assert(DCT);
    }
};
