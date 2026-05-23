#pragma once

#include <cmath>
#include <array>
#include <vector>
#include <fftw3.h>
#include <assert.h>

template <unsigned Dimension>
class Dimensions
{
    template <unsigned Dimension2>
    friend class Function;
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

    unsigned inline points(unsigned i)
    {
        return _points[i];
    }

    size_t inline stride(unsigned i)
    {
        return _stride[i];
    }

    Dimensions(const unsigned *p, bool column_major = false) : column_major(column_major)
    {
        _size = 1;
        if (column_major)
        {
            for (int i = Dimension - 1; i >= 0; i--)
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

    Dimensions(unsigned grid, bool column_major = false) : column_major(column_major)
    {
        _size = 1;
        if (column_major)
        {
            for (int i = Dimension - 1; i >= 0; i--)
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
class Weighty : public Dimensions<Dimension>
{
public:
    // const unsigned dimension = Dimension;
    // const bool column_major;

    typedef std::array<float, Dimension> Event;

    size_t size()
    {
        return Dimensions<Dimension>::size();
    }

    size_t stride(unsigned i)
    {
        return Dimensions<Dimension>::stride(i);
    }

    unsigned points(unsigned i)
    {
        return Dimensions<Dimension>::points(i);
    }

    size_t inline events() const
    {
        return _events;
    }

    class Coordinate : public std::array<unsigned, Dimension>
    {
    private:
        Weighty<Dimension> &weighty;

    public:
        Coordinate inline &operator=(const size_t x)
        {
            for (unsigned i = 0; i < Dimension; i++)
                this[i] = (x / weighty.stride(i)) % weighty.points(i);
            return *this;
        }

        operator size_t() const
        {
            unsigned long x = 0;
            for (unsigned i = 0; i < Dimension; i++)
                x += *this[i] * weighty.stride(i);
            return x;
        }

        Coordinate(Weighty<Dimension> &weighty) : weighty(weighty) {}
    };

    // struct
    // {
    //     template <unsigned Dimension2>
    //     friend class Weighty;

    // private:
    //     unsigned data[Dimension];

    //     unsigned inline &operator[](unsigned i)
    //     {
    //         return data[i];
    //     }

    // public:
    //     unsigned inline operator()(unsigned i) const
    //     {
    //         return data[i];
    //     }
    // } points, stride;

    template <typename Type>
    class Function
    {
        template <unsigned Dimension2>
        friend class Weighty;

        template <unsigned Dimension2>
        friend class Laplacian;

        template <unsigned Dimension2>
        friend class Covariant;

    private:
        Weighty<Dimension> &weighty;
        Type *data;

    protected:
        Type &operator[](size_t x)
        {
            return data[x];
        }

    public:
        Type inline operator()(size_t x)
        {
            return data[x];
        }

        const Type inline *operator()() const
        {
            return data;
        }

        void write(std::ostream &out) const
        {
            out.write(reinterpret_cast<const char *>(Dimensions::_points), Dimension * sizeof(unsigned));
            out.write(reinterpret_cast<const char *>(data), Dimensions::size() * sizeof(Type));
        }

        void write(std::string filename) const
        {
            std::ofstream out(filename, std::ios::binary | std::ios::trunc);
            write(out);
            out.close();
        }

        Function(Weighty<Dimension> &weighty) : weighty(weighty)
        {
            // data = (Type *)fftwf_malloc(weighty.size() * sizeof(Type));
            data = (Type *)malloc(weighty.size() * sizeof(Type));
            Type x = *data;
            *data = x;
        }

        ~Function()
        {
            // fftwf_free(data);
            free(data);
        }
    };

    template <typename Type>
    class FunctionVector : public std::vector<Function<Type>>
    {
    public:
        FunctionVector(Weighty<Dimension> &weighty)
        {
            this->reserve(Dimension);
            for (unsigned i = 0; i < Dimension; i++)
                this->emplace_back(weighty);
        }
    };

    template <typename Type>
    class FunctionMatrix : public std::vector<std::vector<Function<Type>>>
    {
    public:
        FunctionMatrix(Weighty<Dimension> &weighty)
        {
            this->reserve(Dimension);
            for (unsigned i = 0; i < Dimension; i++)
            {
                this->emplace_back(std::vector<Function<Type>>());
                this->back().reserve(Dimension);
                for (unsigned j = 0; j < Dimension; j++)
                    this->back().emplace_back(weighty);
            }
        }
    };

    Function<float> weight = Function<float>(*this);
    Function<float> density = Function<float>(*this);
    Function<float> quantile = Function<float>(*this);
    Function<unsigned short> klass = Function<unsigned short>(*this);

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

    void filter(Function<float> input, Function<float> output, float radius = 0.01f, bool normalize = false)
    {
        // Apply a Gaussian filter to the input function using the DCT.
        // The radius is specified as a fraction of the total size of the grid.
        Function<float> cosine = Function<float>(*this);
        fftwf_execute_r2r((fftwf_plan)DCT, input.data, cosine.data);
        double **kernel = new double *[Dimension];
        for (unsigned i = 0; i < Dimension; i++)
        {
            const double pi = 3.14159265358979323846;
            double *k = kernel[i] = new double[points(i)];
            for (unsigned j = 0; j < points(i); j++)
                k[j] = exp(-2.0 * squared(j * radius * pi));
        }
        for (size_t x = 0; x < size(); x++)
        {
            double k = 1.0;
            for (unsigned i = 0; i < Dimension; i++)
            {
                unsigned j = (x / stride(i)) % points(i);
                k *= kernel[i][j];
            }
            cosine[x] *= k;
        }
        fftwf_execute_r2r((fftwf_plan)DCT, cosine.data, output.data);
        for (unsigned i = 0; i < Dimension; i++)
            free(kernel[i]);
        delete[] kernel;
        if (normalize)
            for (unsigned x = 0; x < size(); x++)
                output[x] /= (float)fft_normalizer;
    }

    void filter(Function<float> data, float radius = 0.01f, bool normalize = false)
    {
        filter(data, data, radius, normalize);
    }

    bool event(const Event &event)
    {
        // Find the bin corresponding to the event
        size_t x = 0;
        double rem[Dimension];
        for (unsigned i = 0; i < Dimension; i++)
        {
            if (event[i] < 0.0f || event[i] >= 1.0f)
                return false;
            rem[i] = event[i] * (points(i) - 1);
            unsigned floor = static_cast<unsigned>(rem[i]);
            rem[i] -= floor;
            x += floor * stride(i);
        }
        // and distribute the weight to the corners of the hypercube containing it.
        for (unsigned neighbor = 0; neighbor < (unsigned)(1 << Dimension); neighbor++)
        {
            size_t offset = 0;
            double weight = 1.0;
            for (unsigned i = 0; i < Dimension; i++)
            {
                if (neighbor & (1 << i))
                {
                    offset += stride(i);
                    weight *= rem[i];
                }
                else
                {
                    weight *= 1.0f - rem[i];
                }
            }
            this->weight[x + offset] += weight;
        }
        ++_events;
        return true;
    }

    unsigned short classify(const Event &event) const
    {
        size_t x = 0;
        for (unsigned i = 0; i < Dimension; i++)
        {
            if (event[i] < 0.0f || event[i] > 1.0f)
                return 0;
            unsigned floor = static_cast<unsigned>(event[i] * (points[i] - 1));
            x += floor * stride[i];
        }
        return (unsigned short)(klass[x]);
    }

    void prepare(float smoothing = .01f)
    {
        if (smoothing > 0.0f)
        {
            filter(weight, density, smoothing);
            float x = density[0];
            density[0] = x;
            for (size_t x = 0; x < size(); x++)
                if (density[x] < 0.0f)
                    density[x] = 0.0f;
        }
        else
        {
            std::copy(weight.data, weight.data + size(), density.data);
        }

        std::vector<float> sorted;
        std::copy(density.data, density.data + size(), std::back_inserter(sorted));
        std::sort(sorted.begin(), sorted.end());
        std::vector<float> summed;
        summed.resize(sorted.size());
        double sum = 0.0;
        for (size_t x = 0; x < size(); x++)
            summed[x] = sum += sorted[x];
        for (size_t x = 0; x < size(); x++)
        {
            quantile[x] = summed.at((std::lower_bound(sorted.begin(), sorted.end(), density[x]) - sorted.begin())) / sum;
            density[x] /= sum;
            // for (unsigned i = 0; i < Dimension; i++)
            // {
            //     unsigned j = (x / stride[i]) % points[i];
            //     // _P[i][j] += density[x];
            // }
        }
    }

    void trim(std::vector<float> &data, float threshold)
    {
        trim(data.data(), threshold);
    }

    void trim(Function<float> func, float threshold)
    {
        trim(func.data, threshold);
    }

    void trim(float *data, float threshold)
    {
        for (size_t x = 0; x < size(); x++)
            if (quantile[x] < threshold)
                data[x] = std::numeric_limits<float>::quiet_NaN();
    }

    Weighty(const unsigned *points, bool column_major = false) : Dimensions<Dimension>(column_major)
    {
        init();
    }

    Weighty(unsigned grid, bool column_major = false) : Dimensions<Dimension>(grid, column_major)
    {
        init();
    }

private:
    size_t _events = 0;
    fftw_r2r_kind kind[Dimension];
    void *DCT;
    unsigned long fft_normalizer = 1;

    inline double squared(double x) { return x * x; };

    void init()
    {
        int fftw_n[Dimension];
        for (unsigned i = 0; i < Dimension; i++)
        {
            kind[i] = FFTW_REDFT00;
            fft_normalizer *= 2 * (points(i) - 1);
            fftw_n[i] = points(i);
        }
        if (Dimensions<Dimension>::column_major)
            std::reverse(fftw_n, fftw_n + Dimension);
        DCT = (void *)fftwf_plan_r2r(Dimension, fftw_n, weight.data, density.data, kind, 0);
        assert(DCT);
    }
};