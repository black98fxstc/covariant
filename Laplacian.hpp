#pragma once

#include <cmath>
#include <array>
#include <vector>
#include "Weighty.hpp"

template <unsigned Dimension>
class Laplacian : public Weighty<Dimension>
{
private:
    class Hypercube
    {
    private:
        Laplacian<Dimension> &laplace;

    public:
        enum State
        {
            unknown,
            outlier,
            ambiguous,
            contiguous,
            assigned
        };

        size_t x;

        Hypercube(Laplacian<Dimension> &laplace, size_t x) : laplace(laplace), x(x)
        {
            laplace.status[x] = unknown;
        };

        float density() const
        {
            return laplace.Weighty<Dimension>::density(x);
        }

        float quantile() const
        {
            return laplace.Weighty<Dimension>::quantile(x);
        }

        float lapacian() const
        {
            return laplace.L[x];
        }

        unsigned short &cluster()
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

    Weighty<Dimension>::template Function<float> L = typename Weighty<Dimension>::Function<float>(*this);
    Weighty<Dimension>::template Function<unsigned short> klass = typename Weighty<Dimension>::Function<unsigned short>(*this);
    Weighty<Dimension>::template Function<typename Hypercube::State> status = typename Weighty<Dimension>::Function<typename Hypercube::State>(*this);

public:
    void analyze(float smoothing = .01f, float threshold = 0.001f)
    {
        Weighty<Dimension>::prepare(smoothing);

        Weighty<Dimension>::for_each_line([this](const Weighty<Dimension>::Line &fiber)
                                          { this->second_derivatives(fiber); });

        Weighty<Dimension>::filter(L, std::pow(Weighty<Dimension>::size(), -1.0f / (float)Dimension), true);

        Weighty<Dimension>::trim(L, threshold);
    }

    double differentialEquation()
    {
        return _differential_error;
    }

    unsigned cluster(float threshold = 0.001f, bool grow = false)
    {
        typename Weighty<Dimension>::Coordinate coord(*this);
        unsigned clusters = 0;
        std::vector<Hypercube> cubes;
        cubes.reserve(Weighty<Dimension>::size());
        for (size_t x = 0; x < Weighty<Dimension>::size(); x++)
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
                        if (coord[i] < Weighty<Dimension>::points(i) - 1 && status[x + Weighty<Dimension>::stride(i)] == Hypercube::unknown)
                            status[x + Weighty<Dimension>::stride(i)] = Hypercube::contiguous;
                        if (coord[i] > 0 && status[x - Weighty<Dimension>::stride(i)] == Hypercube::unknown)
                            status[x - Weighty<Dimension>::stride(i)] = Hypercube::contiguous;
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
                        if (coord[i] < Weighty<Dimension>::points(i) - 1 && status[x + Weighty<Dimension>::stride(i)] == Hypercube::unknown)
                        {
                            status[x + Weighty<Dimension>::stride(i)] = Hypercube::contiguous;
                            klass[x + Weighty<Dimension>::stride(i)] = it->cluster();
                        }
                        if (coord[i] > 0 && status[x - Weighty<Dimension>::stride(i)] == Hypercube::unknown)
                        {
                            status[x - Weighty<Dimension>::stride(i)] = Hypercube::contiguous;
                            klass[x - Weighty<Dimension>::stride(i)] = it->cluster();
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

    Laplacian(unsigned grid, bool column_major = false) : Weighty<Dimension>(grid, column_major) {}

private:
    void second_derivatives(const Weighty<Dimension>::Line &x)
    {
        std::vector<float> S(x.points, 0.0f);
        std::vector<float> T(x.points, 0.0f);

        float max;
        max = 0.0f;
        int m = x.points - 1;
        T[m] = 1.0 / x.delta;
        for (unsigned k = 0; k < x.points; k++)
        {
            if (Weighty<Dimension>::density[x[k]] > max)
            {
                max = Weighty<Dimension>::density[x[k]];
                m = k;
            }
        }

        double tt;
        S[m] = 0.0;
        for (unsigned k = m, j; k < x.points - 1;)
        {
            j = k + 1;
            while (Weighty<Dimension>::density[x[j]] <= 0.0 && j < x.points - 1)
                j++;
            if (Weighty<Dimension>::density[x[j]] <= 0.0)
                tt = 1.0 / x.delta;
            else
                tt = -2.0 * ((std::log(Weighty<Dimension>::density[x[j]]) - std::log(Weighty<Dimension>::density[x[k]])) / squared(x.delta * (j - k)) - S[x[k]] / x.delta / (j - k));
            while (k < j)
            {
                T[k] = (float)tt;
                if (tt > t_max)
                    t_max = tt;
                if (tt < t_min)
                    t_min = tt;
                if (k != x.points - 1)
                {
                    double ss = -tt * x.delta + S[x[k]];
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
            while (j > 0 && Weighty<Dimension>::density[x[j]] <= 0.0)
                j--;
            if (Weighty<Dimension>::density[x[j]] <= 0.0 || Weighty<Dimension>::density[x[k]] <= 0.0)
                tt = 1.0 / squared(x.delta);
            else
                tt = 2.0 * ((std::log(Weighty<Dimension>::density[x[j]]) - std::log(Weighty<Dimension>::density[x[k]])) / squared(x.delta * (k - j)) - S[x[k]] / x.delta / (k - j));
            if (tt > t_max)
                t_max = tt;
            if (tt < t_min)
                t_min = tt;
            while (k > j)
            {
                T[k - 1] = (float)tt;
                double ss = tt * x.delta + S[x[k]];
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

    inline double squared(double x) { return x * x; };

    float s_max = std::numeric_limits<float>::lowest(), s_min = std::numeric_limits<float>::max();
    float t_max = std::numeric_limits<float>::lowest(), t_min = std::numeric_limits<float>::max();
    double _differential_error = 0.0;
};