#pragma once

#include <vector>
#include <array>
#include <string>
#include <fstream>
#include <iostream>
#include <assert.h>
#include <fftw3.h>

// Functions on grids that play nice with FFTW and represent various tensors
template <unsigned Dimension, typename Type>
class Function
{
    template <unsigned Dimension2>
    friend class Weighty;

    template <unsigned Dimension2>
    friend class Laplacian;

    template <unsigned Dimension2>
    friend class Covariant;

private:
    Dimensions<Dimension> &dimensions;
    Type *data;

protected:
    Type inline operator()(size_t x)
    {
        return data[x];
    }

public:
    Type &operator[](size_t x)
    {
        return data[x];
    }

    void write(std::ostream &out) const
    {
        out.write(reinterpret_cast<const char *>(dimensions.points()), Dimension * sizeof(unsigned));
        out.write(reinterpret_cast<const char *>(data), dimensions.size() * sizeof(Type));
    }

    void write(std::string filename) const
    {
        std::ofstream out(filename, std::ios::binary | std::ios::trunc);
        write(out);
        out.close();
    }

    Function(Dimensions<Dimension> &dimensions) : dimensions(dimensions), data(nullptr)
    {
        // The () at the end value-initializes the array (zeroes it for primitives)
        data = new Type[dimensions.size()](); 
        assert(data != nullptr);
    }

    ~Function()
    {
        delete[] data;
    }

    // Rule of Five: Handle move and disable copy to prevent double-free
    Function(const Function&) = delete;
    Function& operator=(const Function&) = delete;

    Function(Function&& other) noexcept : dimensions(other.dimensions), data(other.data) { other.data = nullptr; }
    Function& operator=(Function&&) = delete; // Cannot rebind reference member
};

/**
 * Partial specialization for float types using FFTW's aligned memory allocation.
 */
template <unsigned Dimension>
class Function<Dimension, float>
{
    template <unsigned Dimension2>
    friend class Weighty;

    template <unsigned Dimension2>
    friend class Laplacian;

    template <unsigned Dimension2>
    friend class Covariant;

private:
    Dimensions<Dimension> &dimensions;
    float *data;

protected:
    float inline operator()(size_t x)
    {
        return data[x];
    }

public:
    float &operator[](size_t x)
    {
        return data[x];
    }

    void write(std::ostream &out) const
    {
        out.write(reinterpret_cast<const char *>(dimensions.points()), Dimension * sizeof(unsigned));
        out.write(reinterpret_cast<const char *>(data), dimensions.size() * sizeof(float));
    }

    void write(std::string filename) const
    {
        std::ofstream out(filename, std::ios::binary | std::ios::trunc);
        write(out);
        out.close();
    }

    Function(Dimensions<Dimension> &dimensions) : dimensions(dimensions), data(nullptr)
    {
        data = static_cast<float*>(fftwf_malloc(dimensions.size() * sizeof(float)));
        assert(data != nullptr);
        std::fill(data, data + dimensions.size(), 0.0f);
    }

    ~Function()
    {
        if (data) {
            fftwf_free(data);
        }
    }

    Function(const Function&) = delete;
    Function& operator=(const Function&) = delete;

    Function(Function&& other) noexcept : dimensions(other.dimensions), data(other.data) { other.data = nullptr; }
    Function& operator=(Function&&) = delete;
};

template <unsigned Dimension, typename Type>
class FunctionVector : public std::vector<Function<Dimension, Type>>
{
    template <unsigned Dimension2>
    friend class Weighty;

    template <unsigned Dimension2>
    friend class Laplacian;

    template <unsigned Dimension2>
    friend class Covariant;

public:
    FunctionVector(Dimensions<Dimension> &dimensions)
    {
        this->reserve(Dimension);
        for (unsigned i = 0; i < Dimension; i++)
            this->emplace_back(dimensions);
    }
};

template <unsigned Dimension, typename Type>
class FunctionMatrix : public std::vector<std::vector<Function<Dimension, Type>>>
{
    template <unsigned Dimension2>
    friend class Weighty;

    template <unsigned Dimension2>
    friend class Laplacian;

    template <unsigned Dimension2>
    friend class Covariant;

public:
    FunctionMatrix(Dimensions<Dimension> &dimensions)
    {
        this->reserve(Dimension);
        for (unsigned i = 0; i < Dimension; i++)
        {
            this->emplace_back(std::vector<Function<Dimension, Type>>());
            this->back().reserve(Dimension);
            for (unsigned j = 0; j < Dimension; j++)
                this->back().emplace_back(dimensions);
        }
    }
};

template <unsigned Dimension>
class Weighty;

template <unsigned Dimension, typename Type>
class MarginalFunction : std::array<std::vector<Type>, Dimension>
{
    template <unsigned Dimension2>
    friend class Weighty;

    template <unsigned Dimension2>
    friend class Laplacian;

    template <unsigned Dimension2>
    friend class Covariant;

public:
    MarginalFunction(Dimensions<Dimension> &dimensions)
    {
        for (unsigned i = 0; i < Dimension; i++)
            (*this)[i].resize(dimensions.points(i));
    }
};
