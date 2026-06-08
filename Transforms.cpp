#include "Transforms.hpp"
#include <cmath>
#include <limits>
#include <stdexcept>

const double Transform::LN_10 = std::log(10.0);

static double math_ulp(double d)
{
    if (std::isnan(d))
        return std::numeric_limits<double>::quiet_NaN();
    if (std::isinf(d))
        return std::numeric_limits<double>::infinity();
    if (d == 0.0)
        return std::numeric_limits<double>::denorm_min();
    d = std::abs(d);
    return std::nextafter(d, std::numeric_limits<double>::infinity()) - d;
}

static double adjustA(double W, double M, double A, int bins)
{
    if (bins > 0)
    {
        double zero = (W + A) / (M + A);
        zero = std::rint(zero * bins) / bins;
        return (M * zero - W) / (1 - zero);
    }
    return A;
}

// -----------------------------------------------------------------------------
// Arcsinh
// -----------------------------------------------------------------------------

Arcsinh::Arcsinh(double T, double M, double A, int bins) : T(T), M(M), A(A)
{
    if (T <= 0)
        throw TransformParameterException("T is not positive");
    if (M <= 0)
        throw TransformParameterException("M is not positive");
    if (A < 0)
        throw TransformParameterException("A is negative");
    if (A > M)
        throw TransformParameterException("A is greater than M");

    b = (M + A) * LN_10;
    c = A * LN_10;
    a = T / std::sinh(b - c);
}

Arcsinh::Arcsinh(double T, double M, double A) : Arcsinh(T, M, A, 0) {}

double Arcsinh::scale(double value) const
{
    double x = value / a;
    bool negative = x < 0;
    if (negative)
        x = -x;
    double asinhx = std::log(x + std::sqrt(x * x + 1));
    if (negative)
        return (c - asinhx) / b;
    else
        return (asinhx + c) / b;
}

double Arcsinh::inverse(double scale) const
{
    return a * std::sinh(b * scale - c);
}

void Arcsinh::transform(const std::vector<double> &source, std::vector<double> &destination) const
{
    destination.resize(source.size());
    for (size_t i = 0; i < source.size(); ++i)
    {
        destination[i] = this->scale(source[i]);
    }
}

void Arcsinh::transform(const std::vector<float> &source, std::vector<float> &destination) const
{
    destination.resize(source.size());
    for (size_t i = 0; i < source.size(); ++i)
    {
        destination[i] = static_cast<float>(this->scale(source[i]));
    }
}

double Arcsinh::slope(double scale) const
{
    return a * b * std::cosh(b * scale - c);
}

std::vector<double> Arcsinh::axisLabels() const
{
    double log10x = std::ceil(std::log(T) / LN_10 - M);
    double x = std::exp(LN_10 * log10x);
    int np;
    if (x > T)
    {
        x = T;
        np = 1;
    }
    else
    {
        np = static_cast<int>(std::floor(std::log(T) / LN_10 - log10x)) + 1;
    }
    double B = this->inverse(0);
    int nn;
    if (x > -B)
    {
        nn = 0;
    }
    else if (x == T)
    {
        nn = 1;
    }
    else
    {
        nn = static_cast<int>(std::floor(std::log(-B) / LN_10 - log10x)) + 1;
    }

    std::vector<double> label(nn + np + 1);
    label[nn] = 0;
    for (int i = 1; i <= nn; ++i)
    {
        label[nn - i] = -x;
        label[nn + i] = x;
        x *= 10;
    }
    for (int i = nn + 1; i <= np; ++i)
    {
        label[nn + i] = x;
        x *= 10;
    }
    return label;
}

// -----------------------------------------------------------------------------
// Hyperlog
// -----------------------------------------------------------------------------

Hyperlog::Hyperlog(double T, double W, double M, double A, int bins)
    : T(T), W(W), M(M), A(adjustA(W, M, A, bins))
{

    double actualA = this->A;

    if (this->T <= 0)
        throw TransformParameterException("T is not positive");
    if (this->W < 0)
        throw TransformParameterException("W is negative");
    if (this->M <= 0)
        throw TransformParameterException("M is not positive");
    if (this->W <= 0)
        throw TransformParameterException("W is not positive");
    if (2 * this->W > this->M)
        throw TransformParameterException("W is too large");
    if (-actualA > this->W || actualA + this->W > this->M - this->W)
        throw TransformParameterException("A is too large");

    w = this->W / (this->M + actualA);
    x2 = actualA / (this->M + actualA);
    x1 = x2 + w;
    x0 = x2 + 2 * w;
    b = (this->M + actualA) * LN_10;

    double e2bx0 = std::exp(b * x0);
    double c_a = e2bx0 / w;
    double f_a = std::exp(b * x1) + c_a * x1;
    a = this->T / ((std::exp(b) + c_a) - f_a);
    c = c_a * a;
    f = f_a * a;

    xTaylor = x1 + w / 4;
    double coef = a * std::exp(b * x1);
    taylor.resize(16);
    for (size_t i = 0; i < taylor.size(); ++i)
    {
        coef *= b / (i + 1);
        taylor[i] = coef;
    }
    taylor[0] += c;

    inverse_x0 = inverse(x0);
}

Hyperlog::Hyperlog(double T, double W, double M, double A) : Hyperlog(T, W, M, A, 0) {}
Hyperlog::Hyperlog(double T, double W, double M) : Hyperlog(T, W, M, 0.0) {}
Hyperlog::Hyperlog(double T, double W) : Hyperlog(T, W, Transform::DEFAULT_DECADES, 0.0) {}

double Hyperlog::taylorSeries(double scale) const
{
    double x = scale - x1;
    double sum = taylor.back() * x;
    for (int i = static_cast<int>(taylor.size()) - 2; i >= 0; --i)
    {
        sum = (sum + taylor[i]) * x;
    }
    return sum;
}

double Hyperlog::scale(double value) const
{
    if (value == 0)
        return x1;

    bool negative = value < 0;
    if (negative)
        value = -value;

    double x;
    if (value < inverse_x0)
        x = x1 + value * w / inverse_x0;
    else
        x = std::log(value / a) / b;

    double tolerance = 3 * math_ulp(1.0);
    if (x > 1)
        tolerance = 3 * math_ulp(x);

    for (int i = 0; i < 10; ++i)
    {
        double ae2bx = a * std::exp(b * x);
        double y;
        if (x < xTaylor)
            y = taylorSeries(x) - value;
        else
            y = (ae2bx + c * x) - (f + value);

        double abe2bx = b * ae2bx;
        double dy = abe2bx + c;
        double ddy = b * abe2bx;

        double delta = y / (dy * (1 - y * ddy / (2 * dy * dy)));
        x -= delta;

        if (std::abs(delta) < tolerance)
        {
            if (negative)
                return 2 * x1 - x;
            else
                return x;
        }
    }

    throw std::logic_error("scale() didn't converge");
}

double Hyperlog::inverse(double scale) const
{
    bool negative = scale < x1;
    if (negative)
        scale = 2 * x1 - scale;

    double inv;
    if (scale < xTaylor)
        inv = taylorSeries(scale);
    else
        inv = (a * std::exp(b * scale) + c * scale) - f;

    if (negative)
        return -inv;
    else
        return inv;
}

void Hyperlog::transform(const std::vector<double> &source, std::vector<double> &destination) const
{
    destination.resize(source.size());
    for (size_t i = 0; i < source.size(); ++i)
    {
        destination[i] = this->scale(source[i]);
    }
}

void Hyperlog::transform(const std::vector<float> &source, std::vector<float> &destination) const
{
    destination.resize(source.size());
    for (size_t i = 0; i < source.size(); ++i)
    {
        destination[i] = static_cast<float>(this->scale(source[i]));
    }
}

double Hyperlog::slope(double scale) const
{
    if (scale < x1)
        scale = 2 * x1 - scale;
    return a * b * std::exp(b * scale) + c;
}

std::vector<double> Hyperlog::axisLabels() const
{
    double p = M - 2 * W;
    double log10x = std::ceil(std::log(T) / LN_10 - p);
    double x = std::exp(LN_10 * log10x);
    int np;
    if (x > T)
    {
        x = T;
        np = 1;
    }
    else
    {
        np = static_cast<int>(std::floor(std::log(T) / LN_10 - log10x)) + 1;
    }
    double B = this->inverse(0);
    int nn;
    if (x > -B)
    {
        nn = 0;
    }
    else if (x == T)
    {
        nn = 1;
    }
    else
    {
        nn = static_cast<int>(std::floor(std::log(-B) / LN_10 - log10x)) + 1;
    }

    std::vector<double> label(nn + np + 1);
    label[nn] = 0;
    for (int i = 1; i <= nn; ++i)
    {
        label[nn - i] = -x;
        label[nn + i] = x;
        x *= 10;
    }
    for (int i = nn + 1; i <= np; ++i)
    {
        label[nn + i] = x;
        x *= 10;
    }

    return label;
}

// -----------------------------------------------------------------------------
// Linear
// -----------------------------------------------------------------------------

Linear::Linear(double T, double A, int bins) : T(T), A(A)
{
    if (T <= 0)
        throw TransformParameterException("T is not positive");
    if (T < A)
        throw TransformParameterException("T is less than A");
    if (A < 0)
        throw TransformParameterException("A is negative");

    a = (T + A);
    b = -A;
}

Linear::Linear(double T, double A) : Linear(T, A, 0) {}
Linear::Linear(double T) : Linear(T, 0.0, 0) {}

double Linear::scale(double value) const
{
    return (value - b) / a;
}

double Linear::inverse(double scale) const
{
    return a * scale + b;
}

void Linear::transform(const std::vector<double> &source, std::vector<double> &destination) const
{
    destination.resize(source.size());
    for (size_t i = 0; i < source.size(); ++i)
    {
        destination[i] = this->scale(source[i]);
    }
}

void Linear::transform(const std::vector<float> &source, std::vector<float> &destination) const
{
    destination.resize(source.size());
    for (size_t i = 0; i < source.size(); ++i)
    {
        destination[i] = static_cast<float>(this->scale(source[i]));
    }
}

double Linear::slope(double scale) const
{
    return a;
}

std::vector<double> Linear::axisLabels() const
{
    double r = T - A;
    double log10x = std::floor(std::log(r) / LN_10);
    double x = std::exp(log10x * LN_10);

    double n = r / x;
    if (n <= 1.2)
        x /= 5.0;
    else if (n <= 1.5)
        x /= 4.0;
    else if (n <= 3)
        x /= 2.0;
    else if (n >= 7)
        x *= 2.0;

    int np = static_cast<int>(std::floor(T / x));
    int nn = static_cast<int>(std::floor(A / x));

    std::vector<double> label(nn + np + 1);
    label[nn] = 0;
    for (int i = 1; i <= nn; ++i)
    {
        label[nn - i] = -x * i;
        label[nn + i] = x * i;
    }
    for (int i = nn + 1; i <= np; ++i)
    {
        label[nn + i] = x * i;
    }

    return label;
}

// -----------------------------------------------------------------------------
// Logarithmic
// -----------------------------------------------------------------------------

Logarithmic::Logarithmic(double T, double M, int bins) : T(T), M(M)
{
    if (T <= 0)
        throw TransformParameterException("T is not positive");
    if (M <= 0)
        throw TransformParameterException("M is not positive");

    b = M * LN_10;
    a = T / std::exp(b);
}

Logarithmic::Logarithmic(double T, double M) : Logarithmic(T, M, 0) {}
Logarithmic::Logarithmic(double T) : Logarithmic(T, Transform::DEFAULT_DECADES, 0) {}

double Logarithmic::scale(double value) const
{
    return std::log(value / a) / b;
}

double Logarithmic::inverse(double scale) const
{
    return a * std::exp(b * scale);
}

void Logarithmic::transform(const std::vector<double> &source, std::vector<double> &destination) const
{
    destination.resize(source.size());
    for (size_t i = 0; i < source.size(); ++i)
    {
        destination[i] = this->scale(source[i]);
    }
}

void Logarithmic::transform(const std::vector<float> &source, std::vector<float> &destination) const
{
    destination.resize(source.size());
    for (size_t i = 0; i < source.size(); ++i)
    {
        destination[i] = static_cast<float>(this->scale(source[i]));
    }
}

double Logarithmic::slope(double scale) const
{
    return a * b * std::exp(b * scale);
}

std::vector<double> Logarithmic::axisLabels() const
{
    double log10x = std::ceil(std::log(T) / LN_10 - M);
    double x = std::exp(LN_10 * log10x);
    int n;
    if (x > T)
    {
        x = T;
        n = 0;
    }
    else
    {
        n = static_cast<int>(std::floor(std::log(T) / LN_10 - log10x));
    }

    std::vector<double> label(n + 1);
    for (int i = 0; i <= n; ++i)
    {
        label[i] = x;
        x *= 10;
    }

    return label;
}

// -----------------------------------------------------------------------------
// Logicle
// -----------------------------------------------------------------------------

Logicle::Logicle(double T, double W, double M, double A, int bins)
    : T(T), W(W), M(M), A(adjustA(W, M, A, bins))
{

    double actualA = this->A;

    if (this->T <= 0)
        throw TransformParameterException("T is not positive");
    if (this->W < 0)
        throw TransformParameterException("W is negative");
    if (this->M <= 0)
        throw TransformParameterException("M is not positive");
    if (2 * this->W > this->M)
        throw TransformParameterException("W is too large");
    if (-actualA > this->W || actualA + this->W > this->M - this->W)
        throw TransformParameterException("A is too large");

    w = this->W / (this->M + actualA);
    x2 = actualA / (this->M + actualA);
    x1 = x2 + w;
    x0 = x2 + 2 * w;
    b = (this->M + actualA) * LN_10;
    d = solve(b, w);

    double c_a = std::exp(x0 * (b + d));
    double mf_a = std::exp(b * x1) - c_a / std::exp(d * x1);
    a = this->T / ((std::exp(b) - mf_a) - c_a / std::exp(d));
    c = c_a * a;
    f = -mf_a * a;

    xTaylor = x1 + w / 4;
    double posCoef = a * std::exp(b * x1);
    double negCoef = -c / std::exp(d * x1);

    taylor.resize(16);
    for (size_t i = 0; i < taylor.size(); ++i)
    {
        posCoef *= b / (i + 1.0);
        negCoef *= -d / (i + 1.0);
        taylor[i] = posCoef + negCoef;
    }
    taylor[1] = 0;
}

Logicle::Logicle(double T, double W, double M, double A) : Logicle(T, W, M, A, 0) {}
Logicle::Logicle(double T, double W, double M) : Logicle(T, W, M, 0.0) {}
Logicle::Logicle(double T, double W) : Logicle(T, W, Transform::DEFAULT_DECADES, 0.0) {}

double Logicle::solve(double b, double w)
{
    if (w == 0)
        return b;

    double tolerance = 2 * math_ulp(b);
    double d_lo = 0;
    double d_hi = b;

    double d = (d_lo + d_hi) / 2.0;
    double last_delta = d_hi - d_lo;
    double delta;

    double f_b = -2 * std::log(b) + w * b;
    double f = 2 * std::log(d) + w * d + f_b;
    double last_f = std::numeric_limits<double>::quiet_NaN();

    for (int i = 1; i < 20; ++i)
    {
        double df = 2 / d + w;

        if (((d - d_hi) * df - f) * ((d - d_lo) * df - f) >= 0 || std::abs(1.9 * f) > std::abs(last_delta * df))
        {
            delta = (d_hi - d_lo) / 2.0;
            d = d_lo + delta;
            if (d == d_lo)
                return d;
        }
        else
        {
            delta = f / df;
            double t = d;
            d -= delta;
            if (d == t)
                return d;
        }

        if (std::abs(delta) < tolerance)
            return d;
        last_delta = delta;

        f = 2 * std::log(d) + w * d + f_b;
        if (f == 0 || f == last_f)
            return d;
        last_f = f;

        if (f < 0)
            d_lo = d;
        else
            d_hi = d;
    }

    throw std::logic_error("exceeded maximum iterations in solve()");
}

double Logicle::slope(double scale) const
{
    if (scale < x1)
        scale = 2 * x1 - scale;
    return a * b * std::exp(b * scale) + c * d / std::exp(d * scale);
}

double Logicle::seriesBiexponential(double scale) const
{
    double x = scale - x1;
    double sum = taylor.back() * x;
    for (int i = static_cast<int>(taylor.size()) - 2; i >= 2; --i)
    {
        sum = (sum + taylor[i]) * x;
    }
    return (sum * x + taylor[0]) * x;
}

double Logicle::scale(double value) const
{
    if (value == 0)
        return x1;

    bool negative = value < 0;
    if (negative)
        value = -value;

    double x;
    if (value < f)
    {
        x = x1 + value / taylor[0];
    }
    else
    {
        x = std::log(value / a) / b;
        if (x < x1)
            x = x1;
    }

    double tolerance = 3 * math_ulp(1.0);
    if (x > 1)
        tolerance = 3 * math_ulp(x);

    for (int i = 0; i < 10; ++i)
    {
        double ae2bx = a * std::exp(b * x);
        double ce2mdx = c / std::exp(d * x);
        double y;
        if (x < xTaylor)
        {
            y = seriesBiexponential(x) - value;
        }
        else
        {
            y = (ae2bx + f) - (ce2mdx + value);
        }

        double abe2bx = b * ae2bx;
        double cde2mdx = d * ce2mdx;
        double dy = abe2bx + cde2mdx;
        double ddy = b * abe2bx - d * cde2mdx;

        double delta = y / (dy * (1 - y * ddy / (2.0 * dy * dy)));
        x -= delta;

        if (std::abs(delta) < tolerance)
        {
            if (negative)
                return 2 * x1 - x;
            else
                return x;
        }
    }

    throw std::logic_error("scale() didn't converge");
}

double Logicle::inverse(double scale) const
{
    bool negative = scale < x1;
    if (negative)
        scale = 2 * x1 - scale;

    double inv;
    if (scale < xTaylor)
    {
        inv = seriesBiexponential(scale);
    }
    else
    {
        inv = (a * std::exp(b * scale) + f) - c / std::exp(d * scale);
    }

    if (negative)
        return -inv;
    else
        return inv;
}

void Logicle::transform(const std::vector<double> &source, std::vector<double> &destination) const
{
    destination.resize(source.size());
    for (size_t i = 0; i < source.size(); ++i)
    {
        destination[i] = this->scale(source[i]);
    }
}

void Logicle::transform(const std::vector<float> &source, std::vector<float> &destination) const
{
    destination.resize(source.size());
    for (size_t i = 0; i < source.size(); ++i)
    {
        destination[i] = static_cast<float>(this->scale(source[i]));
    }
}

std::vector<double> Logicle::axisLabels() const
{
    double p = M - 2 * W;
    double log10x = std::ceil(std::log(T) / LN_10 - p);
    double x = std::exp(LN_10 * log10x);
    int np;
    if (x > T)
    {
        x = T;
        np = 1;
    }
    else
    {
        np = static_cast<int>(std::floor(std::log(T) / LN_10 - log10x)) + 1;
    }
    double B = this->inverse(0);
    int nn;
    if (x > -B)
    {
        nn = 0;
    }
    else if (x == T)
    {
        nn = 1;
    }
    else
    {
        nn = static_cast<int>(std::floor(std::log(-B) / LN_10 - log10x)) + 1;
    }

    std::vector<double> label(nn + np + 1);
    label[nn] = 0;
    for (int i = 1; i <= nn; ++i)
    {
        label[nn - i] = -x;
        label[nn + i] = x;
        x *= 10;
    }
    for (int i = nn + 1; i <= np; ++i)
    {
        label[nn + i] = x;
        x *= 10;
    }

    return label;
}