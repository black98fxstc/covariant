#include <random>
#include <fstream>
#include <vector>
#include <string>
#include <algorithm>
#include <iostream>

template <unsigned Dimension, typename Float = float>
class TestData : public std::vector<typename Covariant<Dimension, Float>::Event>
{
    class RandomEvent
    {
    protected:
        inline std::mt19937 &rng()
        {
            thread_local std::random_device rd;
            thread_local std::seed_seq seq{rd(), rd(), rd(), rd()};
            thread_local std::mt19937 gen(seq);
            thread_local std::shuffle_order_engine<std::mt19937, 256> shuffle_gen(gen);
            return gen;
        }
        typename std::uniform_real_distribution<Float>::param_type stddev_param{0.05f, 0.25f};
        typename std::uniform_real_distribution<Float>::param_type mean_param{0.1f, 0.9f};
        typename std::uniform_real_distribution<Float>::param_type fraction_param{0.0f, 1.0f};
        typename std::uniform_real_distribution<Float>::param_type angle_param{0.0f, 3.14159265358979323846f};
        typename std::uniform_real_distribution<Float>::param_type lambda_param{0.25f, 2.5f};
        typename std::normal_distribution<Float>::param_type halfish_param{0.5f, 0.125f};

        typename std::uniform_int_distribution<unsigned> dimension_distribution{0, Dimension - 1};
        typename std::uniform_real_distribution<Float> uniform_distribution{0.0f, 1.0f};
        typename std::normal_distribution<Float> normal_distribution{0.0f, 1.0f};
        typename std::exponential_distribution<Float> exponential_distribution{1.0f};

    public:
        virtual void sample(typename Covariant<Dimension, Float>::Event &event) = 0;
        virtual ~RandomEvent() = default;
    };

public:
    class Normal : public RandomEvent
    {
    public:
        typename Covariant<Dimension>::Event mean;
        Float stddev;

        void sample(typename Covariant<Dimension, Float>::Event &event)
        {
            for (unsigned i = 0; i < Dimension; i++)
            {
                event[i] = RandomEvent::normal_distribution(RandomEvent::rng(), typename std::normal_distribution<Float>::param_type{mean[i], stddev});
            }
        };

        Normal()
        {
            for (unsigned i = 0; i < Dimension; i++)
                mean[i] = RandomEvent::uniform_distribution(RandomEvent::rng(), RandomEvent::mean_param);
            stddev = RandomEvent::uniform_distribution(RandomEvent::rng(), RandomEvent::stddev_param);
        }
    };

    class Snake : public RandomEvent
    {
    public:
        const double pi = 3.14159265358979323846;
        typename Covariant<Dimension>::Event mean;
        unsigned X, Y;
        Float stddev;
        double s, c;

        void sample(typename Covariant<Dimension, Float>::Event &event)
        {
            double delta = RandomEvent::uniform_distribution(RandomEvent::rng(), RandomEvent::fraction_param);
            double delta_x = delta - 0.5;
            double delta_y = std::sin(delta_x * pi * 0.75) / std::sin(pi * 0.375) / 2;
            for (unsigned i = 0; i < Dimension; i++)
            {
                if (i == X)
                    event[i] = RandomEvent::normal_distribution(RandomEvent::rng(), typename std::normal_distribution<Float>::param_type{(Float)(mean[i] + c * delta_x + s * delta_y), stddev});
                else if (i == Y)
                    event[i] = RandomEvent::normal_distribution(RandomEvent::rng(), typename std::normal_distribution<Float>::param_type{(Float)(mean[i] + c * delta_y - s * delta_x), stddev});
                else
                    event[i] = RandomEvent::normal_distribution(RandomEvent::rng(), typename std::normal_distribution<Float>::param_type{mean[i], stddev});
            }
        };

        Snake()
        {
            X = RandomEvent::dimension_distribution(RandomEvent::rng());
            do
            {
                Y = RandomEvent::dimension_distribution(RandomEvent::rng());
            } while (X == Y);
            for (unsigned i = 0; i < Dimension; i++)
            {
                mean[i] = RandomEvent::uniform_distribution(RandomEvent::rng(), RandomEvent::mean_param);
            }
            double theta = RandomEvent::uniform_distribution(RandomEvent::rng(), RandomEvent::angle_param);
            double scale = RandomEvent::normal_distribution(RandomEvent::rng(), RandomEvent::halfish_param);
            s = std::sin(theta) * scale;
            c = std::cos(theta) * scale;   
            stddev = RandomEvent::uniform_distribution(RandomEvent::rng(), RandomEvent::stddev_param);
        }
    };

    class Exponential : public RandomEvent
    {
        unsigned X;
        Float lambda;
        typename Covariant<Dimension>::Event mean;
        Float stddev;

    public:
        void sample(typename Covariant<Dimension, Float>::Event &event)
        {
            for (unsigned i = 0; i < Dimension; i++)
                if (i == X)
                    event[i] = RandomEvent::exponential_distribution(RandomEvent::rng(), typename std::exponential_distribution<Float>::param_type{lambda});
                else
                    event[i] = RandomEvent::normal_distribution(RandomEvent::rng(), typename std::normal_distribution<Float>::param_type{mean[i], stddev});
        };

        Exponential()
        {
            X = RandomEvent::dimension_distribution(RandomEvent::rng());
            lambda = RandomEvent::uniform_distribution(RandomEvent::rng(), RandomEvent::lambda_param);
            for (unsigned i = 0; i < Dimension; i++)
                mean[i] = RandomEvent::uniform_distribution(RandomEvent::rng(), RandomEvent::mean_param);
            stddev = RandomEvent::uniform_distribution(RandomEvent::rng(), RandomEvent::stddev_param);
        }
    };

    class RandomSample : public RandomEvent
    {
        std::vector<RandomEvent *> population;
        std::vector<Float> fractions;

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
                fractions.push_back(RandomEvent::normal_distribution(RandomEvent::rng(), RandomEvent::halfish_param));
                return;
            }

            unsigned n = 0;
            Float max = fractions[0];
            for (unsigned i = 1; i < fractions.size(); i++)
            {
                if (fractions[i] - fractions[i - 1] > max)
                {
                    max = fractions[i] - fractions[i - 1];
                    n = i;
                }
            }
            if (n == 0)
                fractions.insert(fractions.begin(), max * RandomEvent::normal_distribution(RandomEvent::rng(), RandomEvent::halfish_param));
            else
                fractions.insert(fractions.begin() + n, fractions[n - 1] + max * RandomEvent::normal_distribution(RandomEvent::rng(), RandomEvent::halfish_param));
        }

        void sample(typename Covariant<Dimension, Float>::Event &event)
        {
            int p = std::upper_bound(fractions.begin(), fractions.end(), RandomEvent::uniform_distribution(RandomEvent::rng(), RandomEvent::fraction_param)) - fractions.begin();
            population[p]->sample(event);
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
        for (size_t i = 0; i < num_events; i++)
            test_sample.sample((*this)[i]);
    }

    void save(const std::string &filename)
    {
        std::ofstream outfile(filename, std::ios::binary | std::ios::trunc);
        size_t num_events = this->size();
        outfile.write(reinterpret_cast<const char *>(&num_events), sizeof(num_events));
        outfile.write(reinterpret_cast<const char *>(this->data()), num_events * sizeof(typename Covariant<Dimension, Float>::Event));
        outfile.close();
    }

    bool load(const std::string &filename)
    {
        std::ifstream infile(filename, std::ios::binary);
        if (!infile)
            return false;

        size_t num_events;
        infile.read(reinterpret_cast<char *>(&num_events), sizeof(num_events));
        this->resize(num_events);
        infile.read(reinterpret_cast<char *>(this->data()), num_events * sizeof(typename Covariant<Dimension, Float>::Event));
        infile.close();
        return true;
    }
};