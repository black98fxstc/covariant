#pragma once

#include <cmath>
#include <array>
#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <functional>
#include <thread>
#include <fftw3.h>
#include <cctype>
#include <assert.h>

#include "Dimensions.hpp"
#include "Events.hpp"

// Utilities for sampleing multi-dimensional data
template <unsigned Dimension>
class Weighty : public Dimensions<Dimension>, public Events<Dimension>
{
private:
    size_t _events = 0;
    fftw_r2r_kind kind[Dimension]; // Array of FFTW transform kinds for each dimension
    fftwf_plan DCT;                // FFTW plan for float data
    unsigned long fft_normalizer = 1;
    std::array<std::vector<double>, Dimension> kernel;

protected:
public:
    const double pi = 3.14159265358979323846;

    Function<Dimension, float> weight = Function<Dimension, float>(*this);
    Function<Dimension, float> density = Function<Dimension, float>(*this);
    Function<Dimension, float> quantile = Function<Dimension, float>(*this);
    Function<Dimension, unsigned short> cluster_id = Function<Dimension, unsigned short>(*this);
    MarginalFunction<Dimension, float> P = MarginalFunction<Dimension, float>(*this);

    bool visualize = false, verbose = false, antialias = true, verify = true;

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

    size_t events() const
    {
        return Weighty<Dimension>::_events;
    }

    void reset()
    {
        _events = 0;
        weight.zero();
        density.zero();
        quantile.zero();
        cluster_id.zero();
        P.zero();
    }

    void filter(Function<Dimension, float> &input, Function<Dimension, float> &output, float radius = 0.01f, bool normalize = false)
    {
        // Apply a Gaussian filter to the input function using the DCT.
        // The radius is specified as a fraction of full scale.
        // DCT because smoothing the even half-wave means no probability spill across the end points
        Function<Dimension, float> cosine = Function<Dimension, float>(*this);
        fftwf_execute_r2r(DCT, input.data, cosine.data);
        for (unsigned i = 0; i < Dimension; i++)
            for (unsigned j = 0; j < points(i); j++)
                kernel[i][j] = exp(-2.0 * squared(j * radius * pi));
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
        fftwf_execute_r2r(DCT, cosine.data, output.data);
        if (normalize)
            for (unsigned x = 0; x < size(); x++)
                output[x] /= (float)fft_normalizer;
    }

    void filter(Function<Dimension, float> &data, float radius = 0.01f, bool normalize = false)
    {
        filter(data, data, radius, normalize);
    }

    bool event(const Event<Dimension> &event)
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
                    weight *= 1.0 - rem[i];
                }
            }
            this->weight[x + offset] += weight;
        }
        ++_events;
        return true;
    }

    bool locate(const Event<Dimension> &event, Coordinate<Dimension> &coord)
    {
        for (unsigned i = 0; i < Dimension; i++)
            if (event[i] < 0.0f || event[i] >= 1.0f)
                return false;
            else
                coord[i] = static_cast<unsigned>(event[i] * (points(i) - 1));
        return true;
    }

    unsigned short classify(const Event<Dimension> &event)
    {
        Coordinate<Dimension> coord(*this);
        if (!locate(event, coord))
            return 0;
        size_t x = coord;
        return (unsigned short)(cluster_id[x]);
    }

    void prepare(float smoothing = .01f)
    {
        if (smoothing > 0.0f)
        {
            filter(weight, density, smoothing);
            for (size_t x = 0; x < size(); x++)
                if (density[x] < 0.0f)
                    density[x] = 0.0f;
        }
        else
        {
            density = weight;
        }

        std::vector<float> sorted(this->size());
        for (unsigned i = 0; i < this->size(); i++)
            sorted[i] = this->density[i];
        std::sort(sorted.begin(), sorted.end());
        std::vector<float> summed(this->size());
        double sum = 0.0;
        for (size_t x = 0; x < size(); x++)
            summed[x] = sum += sorted[x];
        Coordinate coord(*this);
        for (size_t x = 0; x < size(); x++)
        {
            quantile[x] = summed.at((std::lower_bound(sorted.begin(), sorted.end(), density[x]) - sorted.begin())) / sum;
            density[x] /= sum;
            coord = x;
            for (unsigned i = 0; i < Dimension; i++)
                P[i][coord[i]] += density[x];
        }

        if (this->visualize)
        {
            density.write("density.bin");
            quantile.write("quantile.bin");
            density.write("density.bin");
            quantile.write("quantile.bin");
        }
    }

    void trim(float *data, float threshold)
    {
        for (size_t x = 0; x < size(); x++)
            if (quantile[x] < threshold)
                data[x] = std::numeric_limits<float>::quiet_NaN();
    }

    void trim(Function<Dimension, float> &func, float threshold)
    {
        for (size_t x = 0; x < size(); x++)
            if (quantile[x] < threshold)
                func[x] = std::numeric_limits<float>::quiet_NaN();
    }

    Weighty(const unsigned *points, bool column_major = false) : Dimensions<Dimension>(points, column_major)
    {
        init_fftw();
    }

    Weighty(unsigned grid, bool column_major = false) : Dimensions<Dimension>(grid, column_major)
    {
        init_fftw();
    }

    virtual ~Weighty()
    {
        if (DCT) fftwf_destroy_plan((fftwf_plan)DCT);
    }

private:
    double squared(double x) { return x * x; };

    void init_fftw()
    {
        if (fftwf_init_threads())
            fftwf_plan_with_nthreads(std::max(1u, std::thread::hardware_concurrency()));

        for (unsigned i = 0; i < Dimension; i++)
            kernel[i].resize(points(i));

        int fftw_n[Dimension];
        for (unsigned i = 0; i < Dimension; i++)
        {
            kind[i] = FFTW_REDFT00;
            fft_normalizer *= 2 * (points(i) - 1);
            fftw_n[i] = points(i);
        }
        if (this->column_major) // Use 'this->' for clarity when accessing base class members
            std::reverse(fftw_n, fftw_n + Dimension);

        DCT = fftwf_plan_r2r(Dimension, fftw_n, weight.data, density.data, kind, 0);
        assert(DCT);
    }
};