#pragma once
#include <array>
#include <fstream>
#include <iostream>
#include <limits>
#include <cstddef>
#include <vector>
#include <string>
#include <algorithm>

#include "Dimensions.hpp"

template<unsigned Dimension>
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

template<unsigned Dimension>
class Events : public std::vector<Event<Dimension>>
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
                out.write(reinterpret_cast<const char *>(this->data()), count * sizeof(Event<Dimension>));
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

            Event<Dimension> event;
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
                in.read(reinterpret_cast<char *>(this->data()), count * sizeof(Event<Dimension>));
            return in.good();
        }
    }

    Events() {};
};

