#include <cmath>
#include <array>
#include <vector>
#include <fftw3.h>
#include <assert.h>

#include "Dimensions.hpp"
#include "Function.hpp"
#include "Weighty.hpp"
#include "Laplace.hpp"
#include "TestData.hpp"

template <unsigned Dimension>
class Laplace;

template <unsigned Dimension>
class Riemann : public Laplace<Dimension>
{
private:
    float total_R = 0.0f;
    float s_max = std::numeric_limits<float>::lowest(), s_min = std::numeric_limits<float>::max();
    float t_max = std::numeric_limits<float>::lowest(), t_min = std::numeric_limits<float>::max();
    double differential_error = 0.0;
    double factor_error = 0.0;

    FunctionVector<Dimension, float> f = FunctionVector<Dimension, float>(*this);
    FunctionMatrix<Dimension, float> s = FunctionMatrix<Dimension, float>(*this);
    FunctionMatrix<Dimension, float> t = FunctionMatrix<Dimension, float>(*this);
    FunctionVector<Dimension, float> q = FunctionVector<Dimension, float>(*this);
    FunctionVector<Dimension, float> r = FunctionVector<Dimension, float>(*this);
    FunctionVector<Dimension, float> S = FunctionVector<Dimension, float>(*this);
    FunctionVector<Dimension, float> T = FunctionVector<Dimension, float>(*this);

protected:
    Function<Dimension, float> R = Function<Dimension, float>(*this);
    MarginalFunction<Dimension, float> Q = MarginalFunction<Dimension, float>(*this);

public:
    void analyze(float smoothing = .01f, float threshold = 0.001f)
    {
        total_R = 0.0f;
        differential_error = 0.0;
        factor_error = 0.0;
        s_max = std::numeric_limits<float>::lowest();
        s_min = std::numeric_limits<float>::max();
        t_max = std::numeric_limits<float>::lowest();
        t_min = std::numeric_limits<float>::max();
        f.zero();
        s.zero();
        t.zero();
        q.zero();
        r.zero();
        S.zero();
        T.zero();
        R.zero();
        Q.zero();
        this->P.zero();
        this->L.zero();
        this->cluster_id.zero();
        this->status.zero();

        this->prepare(smoothing);

        // Find the basis functions
        this->for_each_line([this](const Line &fiber)
                            { this->basis_functions(fiber); });
        // verify that they factor the density
        if (this->verify)
            this->verify_factorization();

        // Find the natural parameters
        this->for_each_line([this](const Line &fiber)
                            { this->natural_parameters(fiber); });
        // Evaluate the solutions to the differential equations
        if (this->verify)
        {
            this->bounding_box();
            this->for_each_line([this](const Line &fiber)
                                { this->verify_differential(fiber); });
        }

        // Supress ringing
        if (this->antialias)
        {
            for (unsigned j = 0; j < Dimension; j++)
                for (unsigned i = 0; i < Dimension; i++)
                {
                    this->filter(s[i][j], std::pow(this->size(), -1.0f / (float)Dimension), true);
                    this->filter(t[i][j], std::pow(this->size(), -1.0f / (float)Dimension), true);
                }
        }

        // Lots of arcane summations
        Coordinate<Dimension> coord(*this);
        for (size_t x = 0; x < this->size(); x++)
        {
            coord = x;
            double dual;
            for (unsigned i = 0; i < Dimension; i++)
            {
                dual = 1.0;
                for (unsigned j = 0; j < Dimension; j++)
                {
                    S[j][x] += s[i][j][x];
                    T[j][x] += t[i][j][x];

                    if (j == i)
                        continue;
                    if (f[i][x] > 0.0f)
                        r[i][x] += (t[i][j][x] + t[j][i][x]); // / squared(_f[i][x]);
                    for (unsigned k = 0; k < Dimension; k++)
                    {
                        if (k == i || k == j)
                            continue;
                        if (f[j][x] > 0.0f)
                            q[i][x] += (t[j][k][x] + t[k][j][x]); // / squared(_f[j][x]);
                    }
                    dual *= f[j][x];
                }
                R[x] += r[i][x];

                if (this->quantile[x] < threshold)
                    continue;
                Q[i][coord[i]] += q[i][x] * dual;
            }
            for (unsigned i = 0; i < Dimension; i++)
                this->L[x] += T[i][x];
            R[x] *= 0.5f; // _density[x] / (float)_size;
            if (this->quantile[x] >= threshold)
                total_R += R[x];
        }

        // remove outliers
        this->trim(R, threshold);
        for (unsigned i = 0; i < Dimension; i++)
        {
            this->trim(f[i], threshold);
            this->trim(S[i], threshold);
            this->trim(T[i], threshold);
            for (unsigned j = 0; j < Dimension; j++)
            {
                this->trim(s[i][j], threshold);
                this->trim(t[i][j], threshold);
            }
        }

        // files for MATLAB
        if (this->visualize)
        {
            R.write("R.bin");
            this->L.write("L.bin");
            if (this->verbose)
            {
                f[0].write("f1.bin");
                f[1].write("f2.bin");
                t[0][0].write("t11.bin");
                t[0][1].write("t12.bin");
                t[1][0].write("t21.bin");
                t[1][1].write("t22.bin");
            }
        }
    }

    double factorProbability()
    {
        return factor_error;
    }

    Riemann(const unsigned *points, bool column_major = false) : Laplace<Dimension>(points, column_major) {}

    Riemann(unsigned grid = 256, bool column_major = false) : Laplace<Dimension>(grid, column_major) {}

    ~Riemann() {}

private:
    // factoring the density
    void basis_functions(const Line &x)
    {
        double marginal = (this->density[x[0]] + this->density[x[x.points - 1]]) / 2.0;
        for (unsigned i = 1; i < x.points - 1; i++)
            marginal += this->density[x[i]];
        // keep it real in the denominator
        if (marginal < 1.0 / (double)this->events())
            for (unsigned i = 0; i < x.points; i++)
                f[x.d][x[i]] = 0.0f;
        else
            for (unsigned i = 0; i < x.points; i++)
                f[x.d][x[i]] = this->density[x[i]] / marginal;
    }

    // check that it actually factors the density
    void verify_factorization()
    {
        for (size_t x = 0; x < this->size(); x++)
        {
            if (this->density[x] <= 0.0f)
                continue;
            double product = 1.0;
            for (unsigned i = 0; i < Dimension; i++)
                product *= f[i][x];
            double error = std::abs(product - this->density[x]);
            factor_error += error * this->density[x];
        }
    }

    // basically solving a giant second order partial differential equation on the sample distribution
    void natural_parameters(const Line &x)
    {
        for (unsigned i = 0; i < Dimension; i++)
        {
            float max;
            max = 0.0f;
            int m = x.points - 1;
            t[i][x.d][x[m]] = 1.0 / x.delta;
            for (unsigned k = 0; k < x.points; k++)
            {
                if (f[i][x[k]] > max)
                {
                    max = f[i][x[k]];
                    m = k;
                }
            }

            double tt;
            s[i][x.d][x[m]] = 0.0;
            for (unsigned k = m, j; k < x.points - 1;)
            {
                j = k + 1;
                while (f[i][x[j]] <= 0.0 && j < x.points - 1)
                    j++;
                if (f[i][x[j]] <= 0.0)
                    tt = 1.0 / x.delta;
                else
                    tt = -2.0 * ((std::log((double)f[i][x[j]]) - std::log((double)f[i][x[k]])) / squared(x.delta * (j - k)) - s[i][x.d][x[k]] / (x.delta * (j - k)));
                while (k < j)
                {
                    t[i][x.d][x[k]] = (float)tt;
                    if (k != x.points - 1)
                    {
                        double ss = -tt * x.delta + s[i][x.d][x[k]];
                        s[i][x.d][x[k + 1]] = (float)ss;
                    }
                    k++;
                }
            }
            for (unsigned k = m, j; k > 0;)
            {
                j = k - 1;
                while (j > 0 && f[i][x[j]] <= 0.0)
                    j--;
                if (f[i][x[j]] <= 0.0 || f[i][x[k]] <= 0.0)
                    tt = 1.0 / x.delta;
                else
                    tt = 2.0 * ((std::log((double)f[i][x[k]]) - std::log((double)f[i][x[j]])) / squared(x.delta * (k - j)) - s[i][x.d][x[k]] / (x.delta * (k - j)));
                while (k > j)
                {
                    t[i][x.d][x[k - 1]] = (float)tt;
                    double ss = tt * x.delta + s[i][x.d][x[k]];
                    s[i][x.d][x[k - 1]] = (float)ss;
                    k--;
                }
            }
        }
    }

    // Parameters ranges to scale errors
    void bounding_box()
    {
        for (unsigned i = 0; i < Dimension; i++)
            for (unsigned j = 0; j < Dimension; j++)
                for (unsigned x = 0; x < this->size(); x++)
                {
                    if (t[i][j][x] > t_max)
                        t_max = t[i][j][x];
                    if (t[i][j][x] < t_min)
                        t_min = t[i][j][x];
                    if (s[i][j][x] > s_max)
                        s_max = s[i][j][x];
                    if (s[i][j][x] < s_min)
                        s_min = s[i][j][x];
                }
    }

    // check the math
    void verify_differential(const Line &x)
    {
        double tt, s_diff, t_diff;
        for (unsigned int i = 0; i < Dimension; i++)
        {
            for (unsigned k = 0; k < x.points - 1; k++)
            {
                if (f[i][x[k + 1]] > 0 && f[i][x[k]] > 0)
                {
                    tt = -2.0 * ((std::log(f[i][x[k + 1]]) - std::log(f[i][x[k]])) / squared(x.delta) - s[i][x.d][x[k]] / x.delta);
                    t_diff = std::abs((tt - t[i][x.d][x[k]]) / (t_max - t_min));
                    if (t_diff > differential_error)
                        differential_error = t_diff;
                }
                s_diff = std::abs((s[i][x.d][x[k + 1]] - (-t[i][x.d][x[k]] * x.delta + s[i][x.d][x[k]])) / (s_max - s_min));
                if (s_diff > differential_error)
                    differential_error = s_diff;
            }
        }
    }

    double squared(double x) const { return x * x; };
};
