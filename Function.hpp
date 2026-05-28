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
    // needed for FFTW
    template <unsigned Dimension2>
    friend class Weighty;

   template <unsigned Dimension2>
    friend class Laplace;

   template<unsigned Dimension2>
    friend class Riemann;

private:
    Dimensions<Dimension> &dimensions;
    Type *data;

protected:
    Type &operator[](size_t x)
    {
        return data[x];
    }

public:
    const Type &operator[](size_t x) const
    {
        return data[x];
    }

    Function& operator=(const Function& other)
    {
        if (this != &other) std::copy(other.data, other.data + dimensions.size(), data);
        return *this;
    }

    void zero()
    {
        std::fill(data, data + dimensions.size(), Type{});
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
        data = new Type[dimensions.size()];
        assert(data != nullptr);
    }

    ~Function()
    {
        delete[] data;
    }

    // Rule of Five: Handle move and disable copy to prevent double-free
    Function(const Function &) = delete;
    // Function &operator=(const Function &) = delete;

    Function(Function &&other) noexcept : dimensions(other.dimensions), data(other.data) { other.data = nullptr; }
    Function &operator=(Function &&) = delete; // Cannot rebind reference member
};

/**
 * Partial specialization for float types using FFTW's aligned memory allocation.
 */
template <unsigned Dimension>
class Function<Dimension, float>
{
    // needed for FFTW
    template <unsigned Dimension2>
    friend class Weighty;

   template <unsigned Dimension2>
    friend class Laplace;

   template<unsigned Dimension2>
    friend class Riemann;

private:
    Dimensions<Dimension> &dimensions;
    float *data;

protected:
    float &operator[](size_t x)
    {
        return data[x];
    }

public:
    const float &operator[](size_t x) const
    {
        return data[x];
    }

    Function& operator=(const Function<Dimension, float>& other)
    {
        if (this != &other) std::copy(other.data, other.data + dimensions.size(), data);
        return *this;
    }

    void zero()
    {
        std::fill(data, data + dimensions.size(), 0.0f);
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
        data = static_cast<float *>(fftwf_malloc(dimensions.size() * sizeof(float)));
        assert(data != nullptr);
    }

    ~Function()
    {
        if (data)
            fftwf_free(data);
    }

    Function(const Function &) = delete;
    // Function &operator=(const Function &) = delete;

    Function(Function &&other) noexcept : dimensions(other.dimensions), data(other.data) { other.data = nullptr; }
    Function &operator=(Function &&) = delete;
};

template <unsigned Dimension, typename Type>
class FunctionVector : public std::vector<Function<Dimension, Type>>
{
public:
    void zero()
    {
        for (Function<Dimension, Type>&f : *this)
            f.zero();
    }

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
public:
    void zero()
    {
        for (auto & fv : *this)
            for (auto &f : fv)
                f.zero();
    }

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

template <unsigned Dimension, typename Type>
class MarginalFunction : public std::vector<Function<Dimension, Type>>
{
public:
    void zero()
    {
        for (auto &f : *this)
            f.zero();
    }

    MarginalFunction(Dimensions<Dimension> &dimensions)
    {
        this->reserve(Dimension);
        for (unsigned i = 0; i < Dimension; i++)
            this->emplace_back(dimensions);
    }
};
