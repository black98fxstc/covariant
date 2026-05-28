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

    class Event
    {
    private:
        std::array<float, Dimension> values;

    public:
        float &operator[](size_t i) { return values[i]; }
        const float &operator[](size_t i) const { return values[i]; }
        auto begin() { return values.begin(); }
        auto end() { return values.end(); }
        auto begin() const { return values.begin(); }
        auto end() const { return values.end(); }
        float *data() { return values.data(); }
        const float *data() const { return values.data(); }
    };

    class Events : public std::vector<Event>
    {
    public:
        bool write(const std::string &filename, bool ascii = false) const
        {
            if (ascii)
            {
                std::ofstream out(filename);
                if (!out)
                    return false;

                for (const auto &event : *this)
                {
                    for (unsigned i = 0; i < Dimension; ++i)
                        out << (i == 0 ? "" : " ") << event[i];
                    out << "\n";
                }
                return out.good();
            }
            else
            {
                std::ofstream out(filename, std::ios::binary | std::ios::trunc);
                if (!out)
                    return false;

                size_t count = this->size();
                out.write(reinterpret_cast<const char *>(&count), sizeof(count));
                if (count > 0)
                    out.write(reinterpret_cast<const char *>(this->data()), count * sizeof(Event));
                return out.good();
            }
        }

        bool read(const std::string &filename, bool ascii = false)
        {
            if (ascii)
            {
                std::ifstream in(filename);
                if (!in)
                    return false;

                this->clear();

                // Skip leading whitespace and check if the first line looks like a header
                while (in.peek() != EOF && std::isspace(static_cast<unsigned char>(in.peek())))
                    in.ignore();
                int first = in.peek();
                if (first != EOF && !std::isdigit(static_cast<unsigned char>(first)) && first != '-' && first != '+' && first != '.')
                {
                    // Does not start with a digit or sign; likely a header line.
                    std::string dummy;
                    std::getline(in, dummy);
                }

                Event event;
                while (in >> event[0])
                {
                    for (unsigned i = 1; i < Dimension; ++i)
                    {
                        if (!(in >> event[i]))
                            return false; // Incomplete event record
                    }
                    this->push_back(event);
                }
                return true;
            }
            else
            {
                std::ifstream in(filename, std::ios::binary);
                if (!in)
                    return false;

                size_t count = 0;
                in.read(reinterpret_cast<char *>(&count), sizeof(count));
                if (!in)
                    return false;

                this->resize(count);
                if (count > 0)
                    in.read(reinterpret_cast<char *>(this->data()), count * sizeof(Event));
                return in.good();
            }
        }

        Events() {};
    };
};

template <unsigned Dimension>
class Coordinate
{
private:
    Dimensions<Dimension> *dimensions;
    std::array<unsigned, Dimension> indices;

public:
    unsigned &operator[](size_t i) { return indices[i]; }
    unsigned operator[](size_t i) const { return indices[i]; }
    auto begin() { return indices.begin(); }
    auto end() { return indices.end(); }
    auto begin() const { return indices.begin(); }
    auto end() const { return indices.end(); }

    Coordinate &operator=(const size_t x)
    {
        for (unsigned i = 0; i < Dimension; i++)
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

    Coordinate(Dimensions<Dimension> &dimensions) : dimensions(&dimensions), indices{} {}

    // Coordinates are fixed to their specific dimension instance
    Coordinate(const Coordinate<Dimension> &) = delete;
    Coordinate &operator=(const Coordinate<Dimension> &) = delete;
};
