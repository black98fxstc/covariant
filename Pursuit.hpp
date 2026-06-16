#pragma once

#include <cstdint>
#include <vector>

#include <nlohmann/json.hpp>

#include <Dimensions.hpp>
#include <Boundary.hpp>
#include <Modal.hpp>

using json = nlohmann::json;

typedef uint16_t Measurement;
typedef uint32_t Count;
typedef uint32_t Unique;

class EPP_Parameters
{
public:
    // N = 256 gives points and features a precision of roughly two significant figures

    static const unsigned short N = 256; // resolution of points and boundaries
                                    // optimized when there are lots of small factors

    double W = (double)2 / (double)N; // standard deviation of kernel,
                                        // this is the highest achievable resolution for DBM

    enum Goal
    {                    // which objective function
        best_separation, // lowest edge weight
        best_balance     // edge weight biased towards more even splits
    } goal = best_balance;

    struct KLD // KLD threshold for informative cases
    {
        double Normal1D = .04;     // is the measurement just normal
        double Exponential1D = .2; // is this an exponential tail (CyToF)

        KLD(
            double Normal1D = .04,
            double Exponential1D = .2) noexcept
            : Normal1D(Normal1D), Exponential1D(Exponential1D){};
    };

    KLD kld{.04, .2};

    std::vector<Measurement> censor; // omit measurements from consideration

    // algorithm tweaks

    bool recursive = true; // restart process on the two subsets

    Count finalists = 1; // remember this many of the best candidates

    size_t min_events = 0;    // minimum events to try to split, max sigma squared
    double min_relative = 0; // minimum fraction of total events to try to split

    double tolerance = .01; // default tolerance for polygon simplification

    // implementation details, not intended for general use

    double sigma = 3;        // threshold for starting a new cluste
    Count max_clusters = 12; // most clusters the graph logic should handle

    double kernelRadius(unsigned int pass) const noexcept
    {
        return this->W * pow(sqrt(2.0), pass) / sqrt(2.0);
    }

    bool isCensored(Measurement measurement) const noexcept
    {
        return std::find(this->censor.begin(), this->censor.end(), measurement) != this->censor.end();
    }

    explicit operator json() const noexcept;

    EPP_Parameters &operator=(const json &encoded);

    EPP_Parameters(const json &encoded) { *this = encoded; };

    EPP_Parameters(
        Goal goal = best_balance,
        KLD kld = {.04, .2},
        double W = 2.0 / (double)N)
        : W(W), goal(goal), kld(kld), censor(0), finalists(1),
            sigma(3), max_clusters(12), tolerance(.01){};
};

const EPP_Parameters Default;

enum EPP_Status
{
    EPP_success,
    EPP_characterized,
    EPP_no_cluster,
    EPP_not_significant,
    EPP_threshold,
    EPP_error
};

const static std::vector<std::string> Status_strings{
    "success", "characterized", "no_cluster", "not_significant", "threshold", "error"};

    /**
 * structures defining the result of an ananlysis
 **/
struct Point
{
    Coordinate i, j;

    inline double x() const noexcept { return (double)i / (double)EPP_Parameters::N; };
    inline double y() const noexcept { return (double)j / (double)EPP_Parameters::N; };

    inline bool operator==(const Point &other) const noexcept
    {
        return this->i == other.i && this->j == other.j;
    }

    inline bool operator!=(const Point &other) const noexcept
    {
        return !(*this == other);
    }

    Point(Coordinate i, Coordinate j) noexcept : i(i), j(j){};

    Point() noexcept : i(0), j(0){};

    ~Point() = default;
};

class Polygon : public std::vector<Point>
{
public:
    operator json() const noexcept;
};
