#pragma once

#include <cmath>
#include <array>
#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <functional>
#include <fftw3.h>
#include <assert.h>

template <unsigned Dimension>
class Weighty : public Dimensions<Dimension>
{
private:
    size_t _events = 0;
    fftw_r2r_kind kind[Dimension];
    void *DCT;
    unsigned long fft_normalizer = 1;

public:
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
        return Weighty<Dimension>::_events;
    }

    Function<Dimension, float> weight = Function<Dimension, float>(*this);
    Function<Dimension, float> density = Function<Dimension, float>(*this);
    Function<Dimension, float> quantile = Function<Dimension, float>(*this);
    Function<Dimension, unsigned short> klass = Function<Dimension, unsigned short>(*this);

    MarginalFunction<Dimension, float> P = MarginalFunction<Dimension, float>(*this);
    MarginalFunction<Dimension, float> Q = MarginalFunction<Dimension, float>(*this);

    void filter(Function<Dimension, float> &input, Function<Dimension, float> &output, float radius = 0.01f, bool normalize = false)
    {
        // Apply a Gaussian filter to the input function using the DCT.
        // The radius is specified as a fraction of the total size of the grid.
        Function<Dimension, float> cosine = Function<Dimension, float>(*this);
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

    void filter(Function<Dimension, float> &data, float radius = 0.01f, bool normalize = false)
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

    static constexpr Coordinate<Dimension> OOB = Coordinate<Dimension>(nullptr);

    bool locate(const Event &event, Coordinate<Dimension> &coord)
    {
        for (unsigned i = 0; i < Dimension; i++)
            if (event[i] < 0.0f || event[i] > 1.0f)
                return false;
            else
                coord[i] = static_cast<unsigned>(event[i] * (points(i) - 1));
        return coord;
    }

    unsigned short classify(const Event &event)
    {
        Coordinate<Dimension> coord(*this);
        if (!locate(event, coord))
            return false;
        size_t x = coord;
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
        Coordinate coord(*this);
        for (size_t x = 0; x < size(); x++)
        {
            quantile[x] = summed.at((std::lower_bound(sorted.begin(), sorted.end(), density[x]) - sorted.begin())) / sum;
            density[x] /= sum;
            coord = x;
            for (unsigned i = 0; i < Dimension; i++)
                P[i][coord[i]] += density[x];
        }
        density.write("density.bin");
        quantile.write("quantile.bin");
    }

    void trim(std::vector<float> &data, float threshold)
    {
        trim(data.data(), threshold);
    }

    void trim(Function<Dimension, float> &func, float threshold)
    {
        trim(func.data, threshold);
    }

    void trim(float *data, float threshold)
    {
        for (size_t x = 0; x < size(); x++)
            if (quantile[x] < threshold)
                data[x] = std::numeric_limits<float>::quiet_NaN();
    }

    Weighty(const unsigned *points, bool column_major = false) : Dimensions<Dimension>(points, column_major)
    {
        init();
    }

    Weighty(unsigned grid, bool column_major = false) : Dimensions<Dimension>(grid, column_major)
    {
        init();
    }

private:
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