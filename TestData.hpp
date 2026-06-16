#pragma once

#include <random>
#include <fstream>
#include <vector>
#include <string>
#include <algorithm>
#include <cmath>
#include <iostream>

// #include "Dimensions.hpp"
// #include "Events.hpp"

template <unsigned Dimension>
class Weighty;

class TestData : public std::vector<std::vector<float>>
{
public:
    const unsigned dimension;

    TestData(unsigned dimension) : dimension(dimension) {};

    class RandomEvent
    {
    public:
        const unsigned dimension;

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

        typename std::uniform_int_distribution<unsigned> dimension_distribution{0, dimension - 1};
        typename std::uniform_real_distribution<float> uniform_distribution{0.0f, 1.0f};
        typename std::normal_distribution<float> normal_distribution{0.0f, 1.0f};
        typename std::exponential_distribution<float> exponential_distribution{1.0f};

    public:
        virtual unsigned short sample(std::vector<float> &event) = 0;
        RandomEvent(unsigned dimension) : dimension(dimension) {};
        virtual ~RandomEvent() = default;
    };

    std::vector<unsigned short> _truth;

public:
    class Normal : public RandomEvent
    {
    public:
        std::vector<float> mean;
        float stddev;

        unsigned short sample(std::vector<float> &event) override
        {
            for (unsigned i = 0; i < dimension; i++)
                event[i] = this->normal_distribution(this->rng(), typename std::normal_distribution<float>::param_type{mean[i], stddev});
            return 0;
        };

        Normal(unsigned dimension) : RandomEvent(dimension)
        {
            mean.resize(dimension);
            stddev = this->uniform_distribution(this->rng(), this->stddev_param);
            for (unsigned i = 0; i < dimension; i++)
                mean[i] = this->uniform_distribution(this->rng(), typename std::uniform_real_distribution<float>::param_type{2.0f * stddev, 1.0f - 2.0f * stddev});
        }
    };

    class Snake : public RandomEvent
    {
    public:
        const double pi = 3.14159265358979323846;
        const float half_pi = (float)(pi / 2);
        std::vector<float> mean;
        unsigned X, Y;
        float stddev;
        double s, c;

        unsigned short sample(std::vector<float> &event) override
        {
            // two quarter arcs of the circle, one flipped
            double delta_x, delta_y;
            double theta = this->uniform_distribution(this->rng(), typename std::uniform_real_distribution<float>::param_type{-half_pi, half_pi});
            if (theta > 0)
                delta_x = .5 - .5 * std::cos(theta);
            else
                delta_x = -.5 + .5 * std::cos(theta);
            delta_y = .5 * std::sin(theta);
            for (unsigned i = 0; i < dimension; i++)
                if (i == X)
                    event[i] = this->normal_distribution(this->rng(), typename std::normal_distribution<float>::param_type{(float)(mean[i] + c * delta_x + s * delta_y), stddev});
                else if (i == Y)
                    event[i] = this->normal_distribution(this->rng(), typename std::normal_distribution<float>::param_type{(float)(mean[i] + c * delta_y - s * delta_x), stddev});
                else
                    event[i] = this->normal_distribution(this->rng(), typename std::normal_distribution<float>::param_type{mean[i], stddev});
            return 0;
        };

        Snake(unsigned dimension) : RandomEvent(dimension)
        {
            mean.resize(dimension);
            // random and different dimensions
            X = this->dimension_distribution(this->rng());
            do
            {
                Y = this->dimension_distribution(this->rng());
            } while (X == Y);

            for (unsigned i = 0; i < dimension; i++)
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
        std::vector<float> mean;
        float stddev;

    public:
        unsigned short sample(std::vector<float> &event) override
        {
            for (unsigned i = 0; i < dimension; i++)
                if (i == X)
                    event[i] = this->exponential_distribution(this->rng(), typename std::exponential_distribution<float>::param_type{lambda});
                else
                    event[i] = this->normal_distribution(this->rng(), typename std::normal_distribution<float>::param_type{mean[i], stddev});
            return 0;
        };

        Exponential(unsigned dimension) : RandomEvent(dimension)
        {
            mean.resize(dimension);
            X = this->dimension_distribution(this->rng());
            lambda = this->uniform_distribution(this->rng(), this->lambda_param);
            for (unsigned i = 0; i < dimension; i++)
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

        unsigned short sample(std::vector<float> &event) override
        {
            // pick a random population with appropriate frequency and then sample it
            int p = std::upper_bound(fractions.begin(), fractions.end(), this->uniform_distribution(this->rng(), this->fraction_param)) - fractions.begin();
            population[p]->sample(event);
            return p;
        }

        RandomSample(unsigned dimension) : RandomEvent(dimension) {};

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
        {
            (*this)[i].resize(dimension);
            _truth[i] = test_sample.sample((*this)[i]);
        }
    }

    bool write(const std::string &filename, bool ascii = false) const
    {
        if (ascii)
        {
            std::ofstream out(filename);
            if (!out)
                return false;

            for (const auto &event : *this)
            {
                for (unsigned i = 0; i < dimension; ++i)
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
            {
                for (size_t i = 0; i < count; ++i)
                    out.write(reinterpret_cast<const char *>((*this)[i].data()), dimension * sizeof(float));
            }
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

            std::vector<float> event;
            while (in >> event[0])
            {
                for (unsigned i = 1; i < dimension; ++i)
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
            {
                for (size_t i = 0; i < count; ++i) {
                    (*this)[i].resize(dimension);
                    in.read(reinterpret_cast<char *>((*this)[i].data()), dimension * sizeof(float));
                }
            }
            return in.good();
        }
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