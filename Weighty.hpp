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
#include <mutex>
#include <memory>
#include <cctype>
#include <map>
#include <fftw3.h>
#include <assert.h>
#include <filesystem>
#include <cstdlib>

#include <Eigen/Dense>

#include "Dimensions.hpp"
#include "Function.hpp"
#include "Events.hpp"
#include "Gating.hpp"

// Utilities for sampling multi-dimensional data

const double pi = 3.14159265358979323846;

inline std::mutex& get_fftw_mutex() {
    static std::mutex mtx;
    return mtx;
}

template <unsigned Dimension>
class Weighty : public Dimensions<Dimension>, public Events<Dimension>
{
private:
    size_t _events = 0;
    fftw_r2r_kind kind[Dimension]; // Array of FFTW transform kinds for each dimension
    unsigned long fft_normalizer = 1;
    std::shared_ptr<void> DCT;     // FFTW plan for float data
    inline static std::map<std::pair<std::array<unsigned, Dimension>, bool>, std::shared_ptr<void>> plans;


protected:
public:
    std::unique_ptr<Kernel<Dimension, float>> kernel = std::make_unique<Kernel<Dimension, float>>(*this);
    Function<Dimension, float> weight = Function<Dimension, float>(*this);
    Function<Dimension, float> cosine = Function<Dimension, float>(*this);
    Function<Dimension, float> filtered = Function<Dimension, float>(*this);
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

    void transform (Function<Dimension, float> &input, Function<Dimension, float> &output)
    {
        fftwf_execute_r2r((fftwf_plan)DCT.get(), input.data, output.data);
    }

    void apply_kernel(Function<Dimension, float> &cosine, Function<Dimension, float> &filtered, Kernel<Dimension, float> *kernel)
    {
        Eigen::Map<Eigen::ArrayXf>(filtered.data, size()) = Eigen::Map<Eigen::ArrayXf>(cosine.data, size()) * Eigen::Map<const Eigen::ArrayXf>(kernel->data, size());
    }

    void filter(Function<Dimension, float> &input, Function<Dimension, float> &output, Kernel<Dimension, float> *kernel, bool normalize = false)
    {
        // DCT because smoothing the even half-wave means no probability spill across the end points
        fftwf_execute_r2r((fftwf_plan)DCT.get(), input.data, cosine.data);
        Eigen::Map<Eigen::ArrayXf>(cosine.data, size()) *= Eigen::Map<const Eigen::ArrayXf>(kernel->data, size());
        fftwf_execute_r2r((fftwf_plan)DCT.get(), cosine.data, output.data);
        if (normalize)
        {
            float inv_norm = 1.0f / static_cast<float>(fft_normalizer);
            Eigen::Map<Eigen::ArrayXf>(output.data, size()) *= inv_norm;
        }
    }

    void filter(Function<Dimension, float> &input, Function<Dimension, float> &output, float radius = 0.01f, bool normalize = false)
    {
        kernel->radius(radius);
        filter(input, output, kernel.get(), normalize);
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

    bool locate(const Event<Dimension> &event, Coordinates<Dimension> &coord)
    {
        for (unsigned i = 0; i < Dimension; i++)
            if (event[i] < 0.0f || event[i] >= 1.0f)
                return false;
            else
                coord[i] = static_cast<unsigned>(event[i] * (points(i) - 1));
        return true;
    }

    unsigned short classify(const Coordinates<Dimension> &coord)
    {
        size_t x = coord;
        return cluster_id[x];
    }

    unsigned short classify(const Event<Dimension> &event)
    {
        Coordinates<Dimension> coord(*this);
        if (!locate(event, coord))
            return 0;
        return classify(coord);
    }

    void prepare(float smoothing = .01f)
    {
        if (smoothing > 0.0f)
        {
            filter(weight, density, smoothing, false);
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
        Coordinates coord(*this);
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

    virtual ~Weighty() = default;

private:
    void init_fftw()
    {
        std::lock_guard<std::mutex> lock(get_fftw_mutex());

        int fftw_n[Dimension];
        std::array<unsigned, Dimension> pts;
        for (unsigned i = 0; i < Dimension; i++)
        {
            kind[i] = FFTW_REDFT00;
            fft_normalizer *= 2 * (points(i) - 1);
            fftw_n[i] = points(i);
            pts[i] = points(i);
        }
        if (this->column_major)
            std::reverse(fftw_n, fftw_n + Dimension);

        auto key = std::make_pair(pts, this->column_major);
        auto it = plans.find(key);
        if (it != plans.end())
        {
            DCT = it->second;
            return;
        }

        fftwf_plan raw_plan = fftwf_plan_r2r(Dimension, fftw_n, weight.data, density.data, kind, FFTW_MEASURE);
        assert(raw_plan);
        DCT = std::shared_ptr<void>(raw_plan, [](void *p) {
            if (p) fftwf_destroy_plan((fftwf_plan)p);
        });
        plans.insert({key, DCT});
    }
};