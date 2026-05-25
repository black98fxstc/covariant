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

    size_t inline operator[](unsigned i) const
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

    size_t inline size() const
    {
        return _size;
    }

    unsigned inline points(unsigned i) const
    {
        return _points[i];
    }

    const unsigned inline *points() const
    {
        return _points;
    }

    size_t inline stride(unsigned i) const
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
            for (unsigned i = Dimension; i-- > 0; )
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
            for (unsigned i = Dimension; i-- > 0; )
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

template <unsigned Dimension>
class Coordinate : public std::array<unsigned, Dimension>
{
private:
    Dimensions<Dimension> *dimensions;

public:
    Coordinate inline &operator=(const size_t x)
    {
        for (unsigned i = 0; i < Dimension; i++)
            (*this)[i] = (x / dimensions->stride(i)) % dimensions->points(i);
        return *this;
    }

    operator size_t() const
    {
        unsigned long x = 0;
        for (unsigned i = 0; i < Dimension; i++)
            x += (*this)[i] * dimensions->stride(i);
        return x;
    }

    Coordinate<Dimension> &operator=(const Coordinate<Dimension> &other)
    {
        if (this != &other)
        {
            std::copy(other.begin(), other.end(), this->begin());
        }
        return *this;
    }

    Coordinate(Dimensions<Dimension> &dimensions) : dimensions(&dimensions) {}

    // Allows creation of the OOB sentinel
    constexpr Coordinate(std::nullptr_t) : std::array<unsigned, Dimension>{}, dimensions(nullptr) 
    {
        for (auto &val : *this) val = std::numeric_limits<unsigned>::max();
    }
};

