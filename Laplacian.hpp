#include <cmath>
#include <array>
#include <vector>
#include <algorithm>

#include "Weighty.hpp"


template <unsigned Dimension>
class Weighty;

template <unsigned Dimension>
class Laplacian : public Weighty<Dimension>
{
private:
    float s_max = std::numeric_limits<float>::lowest(), s_min = std::numeric_limits<float>::max();
    float t_max = std::numeric_limits<float>::lowest(), t_min = std::numeric_limits<float>::max();
    double _differential_error = 0.0;

    class Hypercube
    {
    private:
        Laplacian<Dimension> &laplace;

    public:
        size_t x;

        enum State
        {
            unknown,
            outlier,
            ambiguous,
            contiguous,
            assigned
        };

        Hypercube(Laplacian<Dimension> &laplace, size_t x) : laplace(laplace), x(x)
        {
            laplace.status[x] = unknown;
        };

        float density() const
        {
            return laplace.density[x];
        }

        float quantile() const
        {
            return laplace.quantile[x];
        }

        float lapacian() const
        {
            return laplace.L[x];
        }

        float &cluster()
        {
            return laplace.klass[x];
        }

        State &status() const
        {
            return laplace.status[x];
        }

        Hypercube &operator=(const Hypercube &other)
        {
            x = other.x;
            return *this;
        }
    };

    Function<Dimension, float> klass = Function<Dimension, float>(*this);
    Function<Dimension, typename Hypercube::State> status = Function<Dimension, typename Hypercube::State>(*this);

protected:
    Function<Dimension, float> L = Function<Dimension, float>(*this);

public:
    void analyze(float smoothing = .01f, float threshold = 0.001f)
    {
        this->prepare(smoothing);

        this->for_each_line([this](const Line &fiber)
            { this->second_derivatives(fiber); });
        this->for_each_line([this](const Line &fiber)
            { this->differential_error(fiber); });

        this->filter(L, std::pow(this->size(), -1.0f / (float)Dimension), true);
    }

    unsigned cluster(float threshold = 0.001f, bool grow = false)
    {
        Coordinate coord(*this);
        unsigned clusters = 0;
        std::vector<Hypercube> cubes;
        cubes.reserve(this->size());
        for (size_t x = 0; x < this->size(); x++)
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
                    size_t x = it->x;
                    coord = x;
                    for (unsigned i = 0; i < Dimension; i++)
                    {
                        if (coord[i] < this->points(i) - 1 && status[x + this->stride(i)] == Hypercube::unknown)
                            status[x + this->stride(i)] = Hypercube::contiguous;
                        if (coord[i] > 0 && status[x - this->stride(i)] == Hypercube::unknown)
                            status[x - this->stride(i)] = Hypercube::contiguous;
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
                it->cluster() = std::numeric_limits<float>::quiet_NaN();
            }
            for (auto it = outliers; it != cubes.end(); ++it)
            {
                it->status() = Hypercube::outlier;
                it->cluster() = std::numeric_limits<float>::quiet_NaN();
            }
        }
        else
        {
            auto contiguous = precision;
            for (auto clustered = cubes.begin(); clustered != cubes.end();)
            {
                for (auto it = clustered; it != contiguous; ++it)
                {
                    size_t x = it->x;
                    coord = x;
                    for (unsigned i = 0; i < Dimension; i++)
                    {
                        if (coord[i] < this->points(i) - 1 && status[x + this->stride(i)] == Hypercube::unknown)
                        {
                            status[x + this->stride(i)] = Hypercube::contiguous;
                            klass[x + this->stride(i)] = it->cluster();
                        }
                        if (coord[i] > 0 && status[x - this->stride(i)] == Hypercube::unknown)
                        {
                            status[x - this->stride(i)] = Hypercube::contiguous;
                            klass[x - this->stride(i)] = it->cluster();
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
        this->trim(L, threshold);
        L.write("laplacian.bin");
        klass.write("classes.bin");
        return clusters;
    }

    double differentialEquation()
    {
        return _differential_error;
    }

    Laplacian(unsigned grid, bool column_major = false) : Weighty<Dimension>(grid, column_major) {}

private:
    void second_derivatives(const Line &x)
    {
        std::vector<float> S(x.points, 0.0f);
        std::vector<float> T(x.points, 0.0f);

        float max;
        max = 0.0f;
        int m = x.points - 1;
        T[m] = 1.0 / x.delta;
        for (unsigned k = 0; k < x.points; k++)
        {
            if (this->density[x[k]] > max)
            {
                max = this->density[x[k]];
                m = k;
            }
        }

        double tt;
        S[m] = 0.0;
        for (unsigned k = m, j; k < x.points - 1;)
        {
            j = k + 1;
            while (this->density[x[j]] <= 0.0 && j < x.points - 1)
                j++;
            if (this->density[x[j]] <= 0.0)
                tt = 1.0 / x.delta;
            else
                tt = -2.0 * ((std::log(this->density[x[j]]) - std::log(this->density[x[k]])) / squared(x.delta * (j - k)) - S[k] / (x.delta * (j - k)));
            while (k < j)
            {
                T[k] = (float)tt;
                if (tt > t_max)
                    t_max = tt;
                if (tt < t_min)
                    t_min = tt;
                if (k != x.points - 1)
                {
                    double ss = -tt * x.delta + S[k];
                    if (ss > s_max)
                        s_max = ss;
                    if (ss < s_min)
                        s_min = ss;
                    S[k + 1] = (float)ss;
                }
                k++;
            }
        }
        for (unsigned k = m, j; k > 0;)
        {
            j = k - 1;
            while (j > 0 && this->density[x[j]] <= 0.0)
                j--;
            if (this->density[x[j]] <= 0.0 || this->density[x[k]] <= 0.0)
                tt = 1.0 / squared(x.delta);
            else
                tt = 2.0 * ((std::log(this->density[x[j]]) - std::log(this->density[x[k]])) / squared(x.delta * (k - j)) - S[k] / (x.delta * (k - j)));
            if (tt > t_max)
                t_max = tt;
            if (tt < t_min)
                t_min = tt;
            while (k > j)
            {
                T[k - 1] = (float)tt;
                double ss = tt * x.delta + S[k];
                if (ss > s_max)
                    s_max = ss;
                if (ss < s_min)
                    s_min = ss;
                S[k - 1] = (float)ss;
                k--;
            }
        }
        for (unsigned k = 0; k < x.points; k++)
            L[x[k]] += T[k];
    }

    void differential_error(const Line &x)
    {
        double tt, s_diff, t_diff;
        std::vector<float> S(x.points, 0.0f);
        std::vector<float> T(x.points, 0.0f);

        for (unsigned k = 0; k < x.points - 1; k++)
        {
            if (this->density[x[k+1]] > 0 && this->density[x[k]] > 0)
            {
                tt = -2.0 * ((std::log(this->density[x[k+1]]) - std::log(this->density[x[k]])) / squared(x.delta) - S[x[k]] / x.delta);
                t_diff = std::abs(T[k] - tt);
                if (t_diff > t_max - t_min)
                    t_diff = 1.0f;
                else
                    t_diff /= t_max - t_min;
            }
            s_diff = std::abs((S[k + 1] - (-T[k] * x.delta + S[k])) / (s_max - s_min));
            if (s_diff > _differential_error)
                _differential_error = s_diff;
        }
    }

    inline double squared(double x) { return x * x; };

};