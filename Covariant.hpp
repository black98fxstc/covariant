#pragma once

#include <cmath>
#include <array>
#include <vector>
#include <fftw3.h>
#include <assert.h>
#include "Weighty.hpp"
#include "Laplacian.hpp"

template <unsigned Dimension>
class Covariant : public Laplacian<Dimension>
{
private:
    float _tot_R = 0.0f;
    float s_max = std::numeric_limits<float>::lowest(), s_min = std::numeric_limits<float>::max();
    float t_max = std::numeric_limits<float>::lowest(), t_min = std::numeric_limits<float>::max();
    double _differential_error = 0.0;
    double _factor_error = 0.0;
    std::vector<float> _R;
    std::vector<float> _cluster;
    // std::vector<typename Hypercube::state> _status;
    std::array<std::vector<float>, Dimension> _P;
    std::array<std::vector<float>, Dimension> _Q;
    std::array<std::vector<float>, Dimension> _var_Q;
    std::array<float, Dimension> _tot_Q;
    std::array<float, Dimension> _var_tot_Q;

    Weighty<Dimension>::template FunctionVector<float> f = typename Weighty<Dimension>::FunctionVector<float>(*this);
    Weighty<Dimension>::template FunctionMatrix<float> s = typename Weighty<Dimension>::FunctionMatrix<float>(*this);
    Weighty<Dimension>::template FunctionMatrix<float> t = typename Weighty<Dimension>::FunctionMatrix<float>(*this);
    Weighty<Dimension>::template FunctionVector<float> q = typename Weighty<Dimension>::FunctionVector<float>(*this);
    Weighty<Dimension>::template FunctionVector<float> r = typename Weighty<Dimension>::FunctionVector<float>(*this);

    Weighty<Dimension>::template FunctionVector<float> S = typename Weighty<Dimension>::FunctionVector<float>(*this);
    Weighty<Dimension>::template FunctionVector<float> T = typename Weighty<Dimension>::FunctionVector<float>(*this);

    Weighty<Dimension>::template Function<float> L = typename Weighty<Dimension>::Function<float>(*this);

public:
    bool inline event(const Weighty<Dimension>::Event &event)
    {
        return Weighty<Dimension>::event(event);
    }

    void analyse(float smoothing = .01f, float threshold = 0.001f)
    {
        Weighty<Dimension>::prepare(smoothing);

        Weighty<Dimension>::for_each_line([this](const Weighty<Dimension>::Line &fiber)
                                          { this->basis_functions(fiber); });
        Weighty<Dimension>::for_each_line([this](const Weighty<Dimension>::Line &fiber)
                                          { this->natural_parameters(fiber); });

        for (size_t x = 0; x < Weighty<Dimension>::size(); x++)
        {
            if (Weighty<Dimension>::density[x] <= 0.0f)
                continue;
            double product = 1.0;
            for (unsigned i = 0; i < Dimension; i++)
                product *= f[i][x];
            double error = std::abs(product - Weighty<Dimension>::density[x]);
            _factor_error += error * Weighty<Dimension>::density[x];
        }
        Weighty<Dimension>::for_each_line([this](const Weighty<Dimension>::Line &fiber)
                                          { this->differential_error(fiber); });

        if (smoothing > 0.0f)
            for (unsigned j = 0; j < Dimension; j++)
                for (unsigned i = 0; i < Dimension; i++)
                {
                    Weighty<Dimension>::filter(s[i][j], std::pow(Weighty<Dimension>::size(), -1.0f / (float)Dimension), true);
                    Weighty<Dimension>::filter(t[i][j], std::pow(Weighty<Dimension>::size(), -1.0f / (float)Dimension), true);
                }

        typename Weighty<Dimension>::Coordinate coord(*this);
        for (size_t x = 0; x < Weighty<Dimension>::size(); x++)
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
                _R[x] += r[i][x];

                if (Weighty<Dimension>::quantile[x] < threshold)
                    continue;
                _Q[i][coord[i]] += q[i][x] * dual;
                _tot_Q[i] += _Q[i][coord[i]] * _P[i][coord[i]];
            }
            for (unsigned i = 0; i < Dimension; i++)
                L[x] += T[i][x];
            _R[x] *= 0.25f; // _density[x] / (float)_size;
            if (Weighty<Dimension>::quantile[x] >= threshold)
                _tot_R += _R[x];
        }

        Weighty<Dimension>::trim(L, threshold);
        Weighty<Dimension>::trim(_R, threshold);
        for (unsigned i = 0; i < Dimension; i++)
        {
            Weighty<Dimension>::trim(f[i], threshold);
            Weighty<Dimension>::trim(S[i], threshold);
            Weighty<Dimension>::trim(T[i], threshold);
            for (unsigned j = 0; j < Dimension; j++)
            {
                Weighty<Dimension>::trim(s[i][j], threshold);
                Weighty<Dimension>::trim(t[i][j], threshold);
            }
        }
    }

    // unsigned bluster(float threshold = 0.001f)
    // {
    //     for (unsigned x = 0; x < Weighty<Dimension>::size(); x++)
    //         if (std::isnan(L[x]))
    //             _cluster[x] = std::numeric_limits<float>::quiet_NaN();
    //         else if (L[x] <= 0.0f)
    //             _cluster[x] = std::numeric_limits<float>::quiet_NaN();
    //         else
    //             _cluster[x] = L[x];
    //     return 0;
    // }

    double factorProbability()
    {
        return _factor_error;
    }

    double differentialEquation()
    {
        return _differential_error;
    }

    Covariant(const unsigned *points, bool column_major = false) : Laplacian<Dimension>(points, column_major) {}

    Covariant(unsigned grid = 256, bool column_major = false) : Laplacian<Dimension>(grid, column_major) {}

    ~Covariant()
    {
    }

private:
    void basis_functions(const Weighty<Dimension>::Line &x)
    {
        double marginal = (Weighty<Dimension>::density(x[0]) + Weighty<Dimension>::density(x[x.points - 1])) / 2.0;
        for (unsigned i = 1; i < x.points - 1; i++)
        {
            marginal += Weighty<Dimension>::density(x[i]);
        }
        if (marginal < 1.0 / (double)Weighty<Dimension>::events())
            for (unsigned i = 0; i < x.points; i++)
                f[x.d][x[i]] = 0.0f;
        else
            for (unsigned i = 0; i < x.points; i++)
                f[x.d][x[i]] = Weighty<Dimension>::density(x[i]) / marginal;
    }

    void natural_parameters(const Weighty<Dimension>::Line &x)
    {
        for (unsigned i = 0; i < Dimension; i++)
        {
            float max;
            max = 0.0f;
            int m = x.points - 1;
            t[i][x.d][m] = 1.0 / x.delta;
            for (unsigned k = 0; k < x.points; k++)
            {
                if (f[i][x[k]] > max)
                {
                    max = f[i][x[k]];
                    m = k;
                }
            }

            double tt;
            s[i][x.d][m] = 0.0;
            for (unsigned k = m, j; k < x.points - 1;)
            {
                j = k + 1;
                while (f[i][x[j]] <= 0.0 && j < x.points - 1)
                    j++;
                if (f[i][x[j]] <= 0.0)
                    tt = 1.0 / x.delta;
                else
                    tt = -2.0 * ((std::log(f[i][x[j]]) - std::log(f[i][x[k]])) / squared(x.delta * (j - k)) - s[i][x.d][x[k]] / x.delta / (j - k));
                while (k < j)
                {
                    t[i][x.d][k] = (float)tt;
                    if (tt > t_max)
                        t_max = tt;
                    if (tt < t_min)
                        t_min = tt;
                    if (k != x.points - 1)
                    {
                        double ss = -tt * x.delta + s[i][x.d][x[k]];
                        if (ss > s_max)
                            s_max = ss;
                        if (ss < s_min)
                            s_min = ss;
                        s[i][x.d][k + 1] = (float)ss;
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
                    tt = 1.0 / squared(x.delta);
                else
                    tt = 2.0 * ((std::log(f[i][x[k]]) - std::log(f[i][x[j]])) / squared(x.delta * (k - j)) - s[i][x.d][x[k]] / x.delta / (k - j));
                if (tt > t_max)
                    t_max = tt;
                if (tt < t_min)
                    t_min = tt;
                while (k > j)
                {
                    t[i][x.d][k - 1] = (float)tt;
                    double ss = tt * x.delta + s[i][x.d][x[k]];
                    if (ss > s_max)
                        s_max = ss;
                    if (ss < s_min)
                        s_min = ss;
                    s[i][x.d][k - 1] = (float)ss;
                    k--;
                }
            }
        }
    }

    void differential_error(const Weighty<Dimension>::Line &x)
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
                    if (t_diff > _differential_error)
                        _differential_error = t_diff;
                }
                s_diff = std::abs((s[i][x.d][x[k + 1]] - (-t[i][x.d][x[k]] * x.delta + s[i][x.d][x[k]])) / (s_max - s_min));
                if (s_diff > _differential_error)
                    _differential_error = s_diff;
            }
        }
    }

    inline double squared(double x) { return x * x; };

    void init(bool column_major)
    {
        for (unsigned i = 0; i < Dimension; i++)
        {
            _P[i].resize(Weighty<Dimension>::points(i));
            _Q[i].resize(Weighty<Dimension>::points(i));
            _var_Q[i].resize(Weighty<Dimension>::points(i));
        }
        _R.resize(Weighty<Dimension>::size());
        // _cluster.resize(Weighty<Dimension>::size());
        // _status.resize(Weighty<Dimension>::size());
    }
};
