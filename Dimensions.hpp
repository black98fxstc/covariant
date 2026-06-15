#pragma once

#include <array>
#include <fstream>
#include <iostream>
#include <limits>
#include <cstddef>

// Structure constants that must be fixed before any other constructors
struct Line
{
    size_t id;
    size_t base;
    size_t stride;
    float delta;
    unsigned d;
    unsigned points;

    size_t operator[](unsigned i) const
    {
        return base + i * stride;
    }
};

template <unsigned Dimension>
class Dimensions
{

protected:
    unsigned _points[Dimension];
    size_t _stride[Dimension];
    size_t _size;

public:
    const unsigned dimension = Dimension;
    const bool column_major;

    size_t size() const
    {
        return _size;
    }

    unsigned points(unsigned i) const
    {
        return _points[i];
    }

    const unsigned *points() const
    {
        return _points;
    }

    size_t stride(unsigned i) const
    {
        return _stride[i];
    }

    void for_each_line(std::function<void(const Line &)> func)
    {
        Line x;
        for (unsigned i = 0; i < Dimension; ++i)
        {
            x.d = i;
            x.points = points(i);
            x.stride = stride(i);
            x.delta = 1.0 / (double)(points(x.d) - 1);

            for (x.id = 0; x.id < size() / points(x.d); x.id++)
            {
                size_t smaller = x.id % stride(x.d);
                size_t larger = x.id / stride(x.d);
                x.base = larger * points(x.d) * stride(x.d) + smaller;
                func(x);
            }
        }
    }

    // Do things that must complete before Weighty construction
    Dimensions(const unsigned *p, bool column_major = false) : _size(1), column_major(column_major)
    {
        if (column_major)
        {
            for (unsigned i = Dimension; i-- > 0;)
            {
                _stride[i] = _size;
                _points[i] = p[i];
                _size *= _points[i];
            }
        }
        else
        {
            for (unsigned i = 0; i < Dimension; i++)
            {
                _stride[i] = _size;
                _points[i] = p[i];
                _size *= _points[i];
            }
        }
    }

    Dimensions(unsigned grid, bool column_major = false) : _size(1), column_major(column_major)
    {
        if (column_major)
        {
            for (unsigned i = Dimension; i-- > 0;)
            {
                _stride[i] = _size;
                _points[i] = grid + 1;
                _size *= _points[i];
            }
        }
        else
        {
            for (unsigned i = 0; i < Dimension; i++)
            {
                _stride[i] = _size;
                _points[i] = grid + 1;
                _size *= _points[i];
            }
        }
    }
};

typedef int16_t Coordinate;

template <unsigned Dimension>
class Coordinates
{
private:
    Dimensions<Dimension> *dimensions;
    std::array<Coordinate, Dimension> indices;

public:
    Coordinate &operator[](size_t i) { return indices[i]; }
    Coordinate operator[](size_t i) const { return indices[i]; }
    auto begin() { return indices.begin(); }
    auto end() { return indices.end(); }
    auto begin() const { return indices.begin(); }
    auto end() const { return indices.end(); }

    Coordinates &operator=(const size_t x)
    {
        for (Coordinate i = 0; i < Dimension; i++)
            indices[i] = (x / dimensions->stride(i)) % dimensions->points(i);
        return *this;
    }

    operator size_t() const
    {
        size_t x = 0;
        for (unsigned i = 0; i < Dimension; i++)
            x += indices[i] * dimensions->stride(i);
        return x;
    }

    Coordinates &operator++()
    {
        for (Coordinate i = 0; i < Dimension; i++)
        {
            if (indices[i] < dimensions->points(i) - 1)
            {
                indices[i]++;
                return *this;
            }
            indices[i] = 0;
        }
        return *this;
    }


    Coordinates(Dimensions<Dimension> &dimensions) : dimensions(&dimensions), indices{} {}

    // Coordinates are fixed to their specific dimension instance
    Coordinates(const Coordinates<Dimension> &) = delete;
    Coordinates &operator=(const Coordinates<Dimension> &) = delete;
};
