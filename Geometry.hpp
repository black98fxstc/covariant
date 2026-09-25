#pragma once

#include <cstdint>
#include <vector>

#include <nlohmann/json.hpp>

#include "Dimensions.hpp"

using json = nlohmann::json;

using Measurement = uint16_t;
using Count = uint32_t;
using Ordinal = uint32_t;

struct Point
{
    Coordinate i, j;

    inline double x() const noexcept { return static_cast<double>(i) / 256.0; }
    inline double y() const noexcept { return static_cast<double>(j) / 256.0; }

    inline bool operator==(const Point &other) const noexcept
    {
        return this->i == other.i && this->j == other.j;
    }

    inline bool operator!=(const Point &other) const noexcept
    {
        return !(*this == other);
    }

    Point(Coordinate i, Coordinate j) noexcept : i(i), j(j) {}
    Point() noexcept : i(0), j(0) {}
};

class Polygon : public std::vector<Point>
{
public:
    operator json() const noexcept;
    static void close_clockwise(Polygon &polygon) noexcept;
    Polygon simplify(double tolerance) noexcept;

private:
    void simplify(double tolerance, Polygon &simplified, size_t lo, size_t hi) const noexcept;
};
