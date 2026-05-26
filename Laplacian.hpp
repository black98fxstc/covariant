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

    FunctionVector<Dimension, float> S = FunctionVector<Dimension, float>(*this);
    FunctionVector<Dimension, float> T = FunctionVector<Dimension, float>(*this);

protected:
    Function<Dimension, float> L = Function<Dimension, float>(*this);
    Function<Dimension, float> klass = Function<Dimension, float>(*this);
    Function<Dimension, typename Hypercube::State> status = Function<Dimension, typename Hypercube::State>(*this);

public:
    bool verify = true;

    void analyze(float smoothing = .01f, float threshold = 0.001f)
    {
        this->prepare(smoothing);

        // find the second derivatives
        this->for_each_line([this](const Line &fiber)
            { this->second_derivatives(fiber); });
        // check the math
        if (this->verify)
        {
            this->bounding_box();
            this->for_each_line([this](const Line &fiber)
                { this->differential_error(fiber); });
        }

        // Supress ringing
        if (this->antialias)
            this->filter(L, std::pow(this->size(), -1.0f / (float)Dimension), true);

        // Remove outliers
        this->trim(L, threshold);

        // Files for MATLAB
        if (this->visualize)
        {
            L.write("laplacian.bin");
            klass.write("classes.bin"); 
        }
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
        // So that clusters are found in decending order by mode
        std::sort(cubes.begin(), precision, [](const Hypercube &a, const Hypercube &b)
                  { return a.density() > b.density(); });
        for (auto clustered = cubes.begin(); clustered != precision;)
        {
            // start with the absolute mode
            clusters++;
            clustered->status() = Hypercube::contiguous;
            auto contiguous = clustered + 1;
            while (true)
            {
                // claim this cube and mark any unknown neighbors as contiguous
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
                // find the newbs
                contiguous = std::partition(clustered, precision, [](const Hypercube &cube)
                                            { return cube.status() == Hypercube::contiguous; });
                if (contiguous == clustered)
                    break;
            }
        }
        // Clean up the unclaimed cubes
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
            // claim any unclaimed neighboring cubes for existing clusters
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

        // Remove outliers
        this->trim(L, threshold);

        // Files for MATLAB
        if (this->visualize)
        {
            L.write("laplacian.bin");
            klass.write("classes.bin");
        }

        return clusters;
    }

    double differentialEquation()
    {
        return _differential_error;
    }

    Laplacian(unsigned grid, bool column_major = false) : Weighty<Dimension>(grid, column_major) {}

private:

    // Compute second derivatibes of the density
    void second_derivatives(const Line &x)
    {
        float max;
        max = 0.0f;
        int m = x.points - 1;
        T[x.d][x[m]] = 1.0 / x.delta;
        for (unsigned k = 0; k < x.points; k++)
        {
            if (this->density[x[k]] > max)
            {
                max = this->density[x[k]];
                m = k;
            }
        }

        // Arcane summs
        double tt;
        S[x.d][x[m]] = 0.0;
        for (unsigned k = m, j; k < x.points - 1;)
        {
            j = k + 1;
            while (this->density[x[j]] <= 0.0 && j < x.points - 1)
                j++;
            if (this->density[x[j]] <= 0.0)
                tt = 1.0 / x.delta;
            else
                tt = -2.0 * ((std::log(this->density[x[j]]) - std::log(this->density[x[k]])) / squared(x.delta * (j - k)) - S[x.d][x[k]] / (x.delta * (j - k)));
            while (k < j)
            {
                T[x.d][x[k]] = (float)tt;
                if (k != x.points - 1)
                {
                    double ss = -tt * x.delta + S[x.d][x[k]];
                    S[x.d][x[k+1]] = (float)ss;
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
                tt = 2.0 * ((std::log(this->density[x[j]]) - std::log(this->density[x[k]])) / squared(x.delta * (k - j)) - S[x.d][x[k]] / (x.delta * (k - j)));
            while (k > j)
            {
                T[x.d][x[k-1]] = (float)tt;
                double ss = tt * x.delta + S[x.d][x[k]];
                S[x.d][x[k-1]] = (float)ss;
                k--;
            }
        }
        for (unsigned k = 0; k < x.points; k++)
            L[x[k]] += T[x.d][x[k]];
    }

    // check the math
    void bounding_box()
    {
        for (unsigned i = 0; i < Dimension; i++)
            for (unsigned x = 0; x < this->size(); x++)
            {
                if (T[i][x] > t_max)
                    t_max = T[i][x];
                if (T[i][x] < t_min)
                    t_min = T[i][x];
                if (S[i][x] > s_max)
                    s_max = S[i][x];
                if (S[i][x] < s_min)
                    s_min = S[i][x];
            }
    }

    void differential_error(const Line &x)
    {
        double tt, s_diff, t_diff;

        for (unsigned k = 0; k < x.points - 1; k++)
        {
            if (this->density[x[k+1]] > 0 && this->density[x[k]] > 0)
            {
                tt = -2.0 * ((std::log(this->density[x[k+1]]) - std::log(this->density[x[k]])) / squared(x.delta) - S[x.d][x[k]] / x.delta);
                t_diff = std::abs(T[x.d][x[k]] - tt) / (t_max - t_min);
                if (t_diff > _differential_error)
                    _differential_error = t_diff;
            }
            s_diff = std::abs((S[x.d][x[k + 1]] - (-T[x.d][x[k]] * x.delta + S[x.d][x[k]])) / (s_max - s_min));
            if (s_diff > _differential_error)
                _differential_error = s_diff;
        }
    }

    inline double squared(double x) { return x * x; };
};