#pragma once

#include <random>
#include <fstream>
#include <vector>
#include <string>
#include <algorithm>
#include <cmath>
#include <iostream>

#include "Dimensions.hpp"

template <unsigned Dimension>
class Weighty;

template <unsigned Dimension>
class TestData : public Dimensions<Dimension>::Events
{
    class RandomEvent
    {
    protected:
        static std::mt19937 &rng()
        {
            static thread_local std::random_device rd;
            static thread_local std::seed_seq seq{rd(), rd(), rd(), rd()};
            static thread_local std::mt19937 gen(seq);
            return gen;
        }
        typename std::uniform_real_distribution<float>::param_type stddev_param{0.0125f, 0.1f};
        typename std::uniform_real_distribution<float>::param_type mean_param{0.1f, 0.9f};
        typename std::uniform_real_distribution<float>::param_type fraction_param{0.0f, 1.0f};
        typename std::uniform_real_distribution<float>::param_type angle_param{0.0f, 3.14159265358979323846f};
        typename std::uniform_real_distribution<float>::param_type lambda_param{0.25f, 2.5f};
        typename std::normal_distribution<float>::param_type halfish_param{0.5f, 0.125f};

        typename std::uniform_int_distribution<unsigned> dimension_distribution{0, Dimension - 1};
        typename std::uniform_real_distribution<float> uniform_distribution{0.0f, 1.0f};
        typename std::normal_distribution<float> normal_distribution{0.0f, 1.0f};
        typename std::exponential_distribution<float> exponential_distribution{1.0f};

    public:
        virtual unsigned short sample(typename Dimensions<Dimension>::Event &event) = 0;
        virtual ~RandomEvent() = default;
    };

    std::vector<unsigned short> _truth;

public:
    class Normal : public RandomEvent
    {
    public:
        typename Dimensions<Dimension>::Event mean;
        float stddev;

        unsigned short sample(typename Dimensions<Dimension>::Event &event)
        {
            for (unsigned i = 0; i < Dimension; i++)
                event[i] = this->normal_distribution(this->rng(), typename std::normal_distribution<float>::param_type{mean[i], stddev});
            return 0;
        };

        Normal()
        {
            stddev = this->uniform_distribution(this->rng(), this->stddev_param);
            for (unsigned i = 0; i < Dimension; i++)
                mean[i] = this->uniform_distribution(this->rng(), typename std::uniform_real_distribution<float>::param_type{2.0f * stddev, 1.0f - 2.0f * stddev});
        }
    };

    class Snake : public RandomEvent
    {
    public:
        const double pi = 3.14159265358979323846;
        const float half_pi = (float)(pi / 2);
        typename Dimensions<Dimension>::Event mean;
        unsigned X, Y;
        float stddev;
        double s, c;

        unsigned short sample(typename Dimensions<Dimension>::Event &event)
        {
            // two quarter arcs of the circle, one flipped
            double delta_x, delta_y;
            double theta = this->uniform_distribution(this->rng(), typename std::uniform_real_distribution<float>::param_type{-half_pi, half_pi});
            if (theta > 0)
                delta_x = .5 - .5 * std::cos(theta);
            else
                delta_x = -.5 + .5 * std::cos(theta);
            delta_y = .5 * std::sin(theta);
            for (unsigned i = 0; i < Dimension; i++)
                if (i == X)
                    event[i] = this->normal_distribution(this->rng(), typename std::normal_distribution<float>::param_type{(float)(mean[i] + c * delta_x + s * delta_y), stddev});
                else if (i == Y)
                    event[i] = this->normal_distribution(this->rng(), typename std::normal_distribution<float>::param_type{(float)(mean[i] + c * delta_y - s * delta_x), stddev});
                else
                    event[i] = this->normal_distribution(this->rng(), typename std::normal_distribution<float>::param_type{mean[i], stddev});
            return 0;
        };

        Snake()
        {
            // random and different dimensions
            X = this->dimension_distribution(this->rng());
            do
            {
                Y = this->dimension_distribution(this->rng());
            } while (X == Y);

            for (unsigned i = 0; i < Dimension; i++)
                mean[i] = this->uniform_distribution(this->rng(), this->mean_param);
            // random orientation and size
            double theta = this->uniform_distribution(this->rng(), this->angle_param);
            double scale = this->normal_distribution(this->rng(), this->halfish_param);
            s = std::sin(theta) * scale;
            c = std::cos(theta) * scale;

            stddev = this->uniform_distribution(this->rng(), this->stddev_param);
        }
    };

    class Exponential : public RandomEvent
    {
        unsigned X;
        float lambda;
        typename Dimensions<Dimension>::Event mean;
        float stddev;

    public:
        unsigned short sample(typename Dimensions<Dimension>::Event &event)
        {
            for (unsigned i = 0; i < Dimension; i++)
                if (i == X)
                    event[i] = this->exponential_distribution(this->rng(), typename std::exponential_distribution<float>::param_type{lambda});
                else
                    event[i] = this->normal_distribution(this->rng(), typename std::normal_distribution<float>::param_type{mean[i], stddev});
            return 0;
        };

        Exponential()
        {
            X = this->dimension_distribution(this->rng());
            lambda = this->uniform_distribution(this->rng(), this->lambda_param);
            for (unsigned i = 0; i < Dimension; i++)
                mean[i] = this->uniform_distribution(this->rng(), this->mean_param);
            stddev = this->uniform_distribution(this->rng(), this->stddev_param);
        }
    };

    class RandomSample : public RandomEvent
    {
        std::vector<RandomEvent *> population;
        std::vector<float> fractions;

    public:
        void subpopulation(RandomEvent *sub)
        {
            population.push_back(sub);

            if (population.size() == 1)
            {
                return;
            }

            if (population.size() == 2)
            {
                fractions.push_back(this->normal_distribution(this->rng(), this->halfish_param));
                return;
            }

            unsigned n = 0;
            float max = fractions[0];
            for (unsigned i = 1; i < fractions.size(); i++)
            {
                if (fractions[i] - fractions[i - 1] > max)
                {
                    max = fractions[i] - fractions[i - 1];
                    n = i;
                }
            }
            if (n == 0)
                fractions.insert(fractions.begin(), max * this->normal_distribution(this->rng(), this->halfish_param));
            else
                fractions.insert(fractions.begin() + n, fractions[n - 1] + max * this->normal_distribution(this->rng(), this->halfish_param));
        }

        unsigned short sample(typename Dimensions<Dimension>::Event &event)
        {
            // pick a random population with appropriate frequency and then sample it
            int p = std::upper_bound(fractions.begin(), fractions.end(), this->uniform_distribution(this->rng(), this->fraction_param)) - fractions.begin();
            population[p]->sample(event);
            return p;
        }

        virtual ~RandomSample()
        {
            for (auto sub : population)
                delete sub;
        }
    };

    void generate(RandomEvent &test_sample, size_t num_events)
    {
        this->resize(num_events);
        _truth.resize(num_events);
        for (size_t i = 0; i < num_events; i++)
            _truth[i] = test_sample.sample((*this)[i]);
    }

    bool writeTruth(const std::string &filename, bool ascii = false) const
    {
        if (ascii)
        {
            std::ofstream out(filename);
            if (!out)
                return false;
            for (auto c : _truth)
                out << c << "\n";
            return out.good();
        }
        else
        {
            std::ofstream out(filename, std::ios::binary | std::ios::trunc);
            if (!out)
                return false;
            out.write(reinterpret_cast<const char *>(_truth.data()), _truth.size() * sizeof(unsigned short));
            return out.good();
        }
    }

    bool readTruth(const std::string &filename, bool ascii = false)
    {
        size_t count = this->size();
        _truth.resize(count);
        if (ascii)
        {
            std::ifstream in(filename);
            if (!in)
                return false;
            for (size_t i = 0; i < count; i++)
            {
                if (!(in >> _truth[i]))
                    return false;
            }
            return true;
        }
        std::ifstream in(filename, std::ios::binary);
        if (!in)
            return false;
        in.read(reinterpret_cast<char *>(_truth.data()), count * sizeof(unsigned short));
        return in.good();
    }

    const std::vector<unsigned short> &truth() const
    {
        return _truth;
    }
};