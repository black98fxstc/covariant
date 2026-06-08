#pragma once

#include <vector>
#include <stdexcept>
#include <string>

class TransformParameterException : public std::invalid_argument
{
public:
    explicit TransformParameterException(const std::string &message)
        : std::invalid_argument(message) {}
};

class Transform
{
public:
    static constexpr double DEFAULT_DECADES = 4.5;
    static const double LN_10;

    virtual ~Transform() = default;

    virtual double scale(double value) const = 0;
    virtual double inverse(double scale) const = 0;
    virtual std::vector<double> axisLabels() const = 0;
    virtual void transform(const std::vector<double> &source, std::vector<double> &destination) const = 0;
    virtual void transform(const std::vector<float> &source, std::vector<float> &destination) const = 0;

protected:
    virtual double slope(double scale) const = 0;
};

class Arcsinh : public Transform
{
public:
    const double T, M, A;
    double a, b, c;

    Arcsinh(double T, double M, double A, int bins);
    Arcsinh(double T, double M, double A);

    double scale(double value) const override;
    double inverse(double scale) const override;
    std::vector<double> axisLabels() const override;
    void transform(const std::vector<double> &source, std::vector<double> &destination) const override;
    void transform(const std::vector<float> &source, std::vector<float> &destination) const override;

protected:
    double slope(double scale) const override;
};

class Hyperlog : public Transform
{
public:
    const double T, W, M, A;
    double a, b, c, f, w, x0, x1, x2;

protected:
    double xTaylor;
    std::vector<double> taylor;
    double inverse_x0;

    Hyperlog(double T, double W, double M, double A, int bins);

public:
    Hyperlog(double T, double W, double M, double A);
    Hyperlog(double T, double W, double M);
    Hyperlog(double T, double W);

    double scale(double value) const override;
    double inverse(double scale) const override;
    std::vector<double> axisLabels() const override;
    void transform(const std::vector<double> &source, std::vector<double> &destination) const override;
    void transform(const std::vector<float> &source, std::vector<float> &destination) const override;

protected:
    double taylorSeries(double scale) const;
    double slope(double scale) const override;
};

class Linear : public Transform
{
public:
    const double T, A;
    double a, b;

protected:
    Linear(double T, double A, int bins);

public:
    Linear(double T, double A);
    explicit Linear(double T);

    double scale(double value) const override;
    double inverse(double scale) const override;
    std::vector<double> axisLabels() const override;
    void transform(const std::vector<double> &source, std::vector<double> &destination) const override;
    void transform(const std::vector<float> &source, std::vector<float> &destination) const override;

protected:
    double slope(double scale) const override;
};

class Logarithmic : public Transform
{
public:
    const double T, M;
    double a, b;

protected:
    Logarithmic(double T, double M, int bins);

public:
    Logarithmic(double T, double M);
    explicit Logarithmic(double T);

    double scale(double value) const override;
    double inverse(double scale) const override;
    std::vector<double> axisLabels() const override;
    void transform(const std::vector<double> &source, std::vector<double> &destination) const override;
    void transform(const std::vector<float> &source, std::vector<float> &destination) const override;

protected:
    double slope(double scale) const override;
};

class Logicle : public Transform
{
public:
    const double T, W, M, A;
    double a, b, c, d, f, w, x0, x1, x2;

protected:
    double xTaylor;
    std::vector<double> taylor;

    Logicle(double T, double W, double M, double A, int bins);
    static double solve(double b, double w);

public:
    Logicle(double T, double W, double M, double A);
    Logicle(double T, double W, double M);
    Logicle(double T, double W);

    double scale(double value) const override;
    double inverse(double scale) const override;
    std::vector<double> axisLabels() const override;
    void transform(const std::vector<double> &source, std::vector<double> &destination) const override;
    void transform(const std::vector<float> &source, std::vector<float> &destination) const override;

protected:
    double slope(double scale) const override;
    double seriesBiexponential(double scale) const;
};