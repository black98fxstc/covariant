#pragma once

#include <vector>
#include <string>
#include <memory>
#include <future>
#include <atomic>
#include <limits>
#include <iostream>
#include <cstdint>

#include <matplot/matplot.h>
#include <nlohmann/json.hpp>

// #include "Pursuit.hpp"
#include "Leonard.hpp"
#include "FlowJo.hpp"
#include "Workers.hpp"
#include "Covariant.hpp"
#include "Weighty.hpp"
#include "Events.hpp"
#include "Gating.hpp"
#include "Samples.hpp"

using json = nlohmann::json;

struct Params
{
    std::vector<std::string> files;
    std::vector<std::string> variables;
    std::vector<std::string> populations;
    float smoothing = 0.01f;
    float threshold = 0.001f;
    unsigned max_clusters = 12;
    size_t min_events = 0;
    float min_relative = 0.0f;
    float kld_normal = 0.04f;
    float kld_exponential = 0.2f;
    float tolerance = 0.01f;
    unsigned grid_size = 256;
    bool antialias = true;
    bool verify = true;
    int analysis_choice = 1;
};

typedef uint16_t Measurement;
typedef uint32_t Count;
typedef uint32_t Ordinal;

struct Point
{
    Coordinate i, j;

    inline double x() const noexcept { return (double)i / (double)256; };
    inline double y() const noexcept { return (double)j / (double)256; };

    inline bool operator==(const Point &other) const noexcept
    {
        return this->i == other.i && this->j == other.j;
    }

    inline bool operator!=(const Point &other) const noexcept
    {
        return !(*this == other);
    }

    Point(Coordinate i, Coordinate j) noexcept : i(i), j(j){};

    Point() noexcept : i(0), j(0){};

    ~Point() = default;
};

class Polygon : public std::vector<Point>
{
public:
    operator json() const noexcept;
};

class Qualify_Results
{
public:
    const Measurement X;
    double KLDn = 0;
    double KLDe = 0;
    bool qualified = false;

    Qualify_Results(const Measurement X) noexcept : X(X) {};
};

class EPP_Node_Results
{
public:
    std::vector<std::vector<float>> quant_data;

    EPP_Node_Results(Weighty<2> &parent) : quant_data(parent.points(0), std::vector<float>(parent.points(1))) {};
};

class Projection_Results
{
public:
    enum Status
    {
        EPP_success,
        EPP_characterized,
        EPP_no_cluster,
        EPP_not_significant,
        EPP_threshold,
        EPP_error
    } outcome = EPP_error;
    const Measurement X, Y;
    Polygon separatrix;
    double score = std::numeric_limits<double>::infinity();
    struct {
        size_t count = 0;
        std::vector<bool> set;
        Polygon polygon;
    } in, out;
    unsigned int pass = 0, clusters = 0, graphs = 0, merges = 0, splits = 0;

    Projection_Results(Measurement X, Measurement Y) noexcept : X(X < Y ? X : Y), Y(X < Y ? Y : X) {};
};

class Pursuit_Results
{
public:
    std::vector<Measurement> qualified;
    double best_score = 0;
    std::unique_ptr<Projection_Results> best_split;
    std::vector<std::future<Pursuit_Results>> future_children;
    std::vector<Pursuit_Results> children;
    std::future<EPP_Node_Results> future_node;
    std::unique_ptr<EPP_Node_Results> EPP_node;

    void wait() noexcept;
};

template <unsigned Dimension>
void for_each_plane(std::function<void(const unsigned i, const unsigned j)> func)
{
    for (unsigned i = 0; i < Dimension - 1; i++)
        for (unsigned j = i + 1; j < Dimension; j++)
            func(i, j);
};

class Laplace_Results;

class Marginal_Results
{
public:
    std::weak_ptr<Laplace_Results> laplace;
    std::vector<std::vector<std::vector<double>>> class_data;
    std::vector<std::vector<double>> quant_data;
    std::string x_label, y_label;
    unsigned clusters;

    Marginal_Results(std::shared_ptr<Laplace_Results> laplace) noexcept : laplace(laplace) {};
};

class Laplace_Results
{
public:
    std::string parent_name;
    unsigned clusters_found;
    unsigned valid_clusters;
    std::vector<std::vector<unsigned>> cluster_events;
    std::ofstream xml_out;
    json report;
    std::vector<size_t> idx;
    std::vector<unsigned short> classification;
    std::vector<std::vector<float>> means;
    std::vector<std::vector<std::vector<float>>> covariances;
    std::vector<std::future<Marginal_Results>> future_marginals;
    std::vector<Marginal_Results> marginals;

    void wait() noexcept
    {
        for (auto &f : future_marginals)
            marginals.push_back(f.get());
    }
};

class Leonard
{
public:
    Params params;
    SelectionState selections;
    size_t laplacian_offset = 0;
    std::atomic<unsigned> laplace_id_counter = 0;
    std::atomic<unsigned> EPP_id_counter = 0;
    Workspace ws;
    std::vector<SampleData> dummy_samples;
    std::vector<std::vector<double>> colors;
    ThreadPool compute_plane{1};//std::thread::hardware_concurrency()};
    ThreadPool control_plane{4};

    int parse_args(int argc, char *argv[]);
    /**
     * @brief Simple HSV to RGB conversion helper.
     * 
     * @param h Hue [0, 1]
     * @param s Saturation [0, 1]
     * @param v Value [0, 1]
     * @return std::vector<double> RGB components
     */

     static std::vector<double> hsv_to_rgb(double h, double s, double v) {
        double r = 0, g = 0, b = 0;
        if (s == 0) {
            r = g = b = v;
        } else {
            double h_pos = (h == 1.0) ? 0.0 : h * 6.0;
            int i = static_cast<int>(std::floor(h_pos));
            double f = h_pos - i;
            double p = v * (1.0 - s);
            double q = v * (1.0 - (s * f));
            double t = v * (1.0 - (s * (1.0 - f)));
            switch (i) {
                case 0: r = v; g = t; b = p; break;
                case 1: r = q; g = v; b = p; break;
                case 2: r = p; g = v; b = t; break;
                case 3: r = p; g = q; b = v; break;
                case 4: r = t; g = p; b = v; break;
                default: r = v; g = p; b = q; break;
            }
        }
        return {r, g, b};
    }

    Qualify_Results do_Qualify(const std::vector<float> *data, const Measurement X, const std::vector<bool> &included, std::string pop_name);

    Projection_Results do_Projection(const std::vector<std::vector<float> *> &data, const Measurement X, const Measurement Y, const std::vector<bool> &included, std::string pop_name);

    EPP_Node_Results do_EPP_Node(const std::vector<std::vector<float>*> data, const std::vector<bool> &included, std::string pop_name);

    Pursuit_Results do_Pursuit(const std::vector<std::vector<float> *> &data, const std::vector<bool> &included, std::string pop_name);

    Marginal_Results do_Marginal(std::shared_ptr<Laplace_Results> laplace, const std::vector<std::vector<float> *> &data, const Measurement i, const Measurement j, std::string pop_name) noexcept
    {
        Marginal_Results results(laplace);

        Weighty<2> marginal(256);
        marginal.visualize = false; // We will save files manually into the dir
        marginal.antialias = params.antialias;
        marginal.verify = params.verify;

        Coordinates<2> marginal_coord(marginal);
        Event<2> marginal_event;
        std::vector<unsigned short> marginal_klass(marginal.points(0) * marginal.points(1), 0);
        std::fill(marginal_klass.begin(), marginal_klass.end(), 0);

        for (unsigned short c = 0; c <= laplace->valid_clusters; ++c)
        {
            for (unsigned e : laplace->cluster_events[c])
            {
                marginal_event[0] = (*data[i])[laplace->idx[e]];
                marginal_event[1] = (*data[j])[laplace->idx[e]];
                marginal.event(marginal_event);
                if (c == 0) continue;

                size_t idx = (size_t)marginal_coord;
                unsigned short d = marginal_klass[idx];
                if (d == 0) marginal_klass[idx] = c;
                else marginal_klass[idx] = std::min(c, d);
            }
        }
        marginal.prepare(params.smoothing);

        using namespace matplot;

        std::vector<std::vector<std::vector<double>>> class_data(3, std::vector<std::vector<double>>(marginal.points(0), std::vector<double>(marginal.points(1))));
        std::vector<std::vector<double>> quant_data(marginal.points(0), std::vector<double>(marginal.points(1)));
        for (unsigned y = 0; y < marginal.points(1); ++y) {
            for (unsigned x = 0; x < marginal.points(0); ++x) {
                size_t idx = x + y * marginal.points(0);
                for (unsigned i = 0; i < 3; ++i)
                    if (marginal_klass[idx] == 0)
                        class_data[i][y][x] = 255;
                    else
                    {
                        unsigned hue = (255 * marginal_klass[idx] / (laplace->valid_clusters + 1));
                        class_data[i][y][x] = 255 * colors[hue][i];
                    }
                quant_data[y][x] = (double)static_cast<const Function<2, float>&>(marginal.quantile)[idx];
            }
        }

        // dispatch plot

        for (unsigned c = 0; c <= laplace->valid_clusters; ++c)
        {
            std::vector<unsigned short> marginal_klass(marginal.size(), 0);

            marginal.reset();
            for (auto &e : laplace->cluster_events[c])
            {
                marginal_event[0] = (*data[i])[laplace->idx[e]];
                marginal_event[1] = (*data[j])[laplace->idx[e]];
                marginal.event(marginal_event);
                if (c == 0) continue;

                size_t idx = (size_t)marginal_coord;
                marginal_klass[idx] = c;
            }
            marginal.prepare(selections.smoothing);

            // make cluster marginal plot
            std::vector<std::vector<std::vector<double>>> class_data(3, std::vector<std::vector<double>>(marginal.points(0), std::vector<double>(marginal.points(1))));
            std::vector<std::vector<double>> quant_data(marginal.points(0), std::vector<double>(marginal.points(1)));

            for (unsigned y = 0; y < marginal.points(1); ++y) {
                for (unsigned x = 0; x < marginal.points(0); ++x) {
                    size_t idx = x + y * marginal.points(1);
                    for (unsigned i = 0; i < 3; ++i)
                        if (marginal_klass[idx] == 0)
                            class_data[i][y][x] = 255;
                        else
                        {
                            unsigned hue = (255 * marginal_klass[idx] / (laplace->valid_clusters + 1));
                            class_data[i][y][x] = 255 * colors[hue][i];
                        }
                    quant_data[y][x] = (double)static_cast<const Function<2, float>&>(marginal.quantile)[idx];
                }
            }

            // dispatch plot

            // auto class_data = std::make_tuple(std::move(class_r), std::move(class_g), std::move(class_b));
            // std::string path = params.img_dir + "/cluster_" + std::to_string(c) + "_" + darwin.labels[i] + "_" + darwin.labels[j] + ".png";
            // darwin.dispatch_plot(path, class_data, quant_data);
            // cluster_images.push_back("images/cluster_" + std::to_string(c) + "_" + darwin.labels[i] + "_" + darwin.labels[j] + ".png");
        }
        // std::string path = params.img_dir + "/sample_" + darwin.labels[i] + "_" + darwin.labels[j] + ".png";
        // darwin.dispatch_plot(path, class_data, quant_data);
        // sample_images.push_back("images/sample_" + darwin.labels[i] + "_" + darwin.labels[j] + ".png");

        return results;
    }

    template <unsigned Dimension>
    std::shared_ptr<Laplace_Results> do_Laplace(const std::vector<std::vector<float> *> &data, const std::vector<bool> &included, std::string pop_name)
    {
        auto results = std::make_shared<Laplace_Results>();
        results->parent_name = pop_name;
        size_t included_events = std::count(included.begin(), included.end(), true);
        if (!included_events)
            return results;

        results->idx.reserve(included_events);
        for (size_t i = 0; i < included.size(); ++i)
            if (included[i])
                results->idx.push_back(i);

        std::cout << "Begin Laplacian clustering on " << results->idx.size() << " events..." << std::endl;
        Events<Dimension> events;
        events.resize(results->idx.size());
        for (unsigned d = 0; d < Dimension; ++d)
            for (size_t i = 0; i < events.size(); ++i)
                events[i][d] = (*data[d])[results->idx[i]];

        Laplace<Dimension> laplace(selections.grid_size);
        laplace.verify = params.verify;
        laplace.antialias = params.antialias;

        size_t valid_events = 0;
        for (const auto &e : events)
            if (laplace.event(e))
                valid_events++;
        std::cout << "Found " << valid_events << " valid events..." << std::endl;

        std::cout << "Calculating the Laplacian of the sample..." << std::endl;
        laplace.analyze(selections.smoothing, selections.threshold);
        if (laplace.differentialError() > .0001)
            std::cout << "Differential equation solution is unusually bad " << laplace.differentialError() << std::endl;
        else
            std::cout << "Consistency checks passed..." << std::endl;

        std::cout << "Performing Laplacian clustering..." << std::endl;
        results->clusters_found = laplace.cluster(selections.threshold);
        results->valid_clusters = std::min(selections.max_clusters, results->clusters_found);
        results->cluster_events.resize(results->valid_clusters + 2);
        results->means.resize(results->valid_clusters + 1);
        results->covariances.resize(results->valid_clusters + 1);
        for (unsigned c = 0; c <= results->valid_clusters; ++c)
        {
            results->means[c].resize(Dimension);
            results->covariances[c].resize(Dimension);
            for (unsigned i = 0; i < Dimension; i++)
                    results->covariances[c][i].resize(Dimension);
        }
        
        Coordinates coord(laplace);
        for (size_t i = 0; i < results->idx.size(); ++i)
        {
            unsigned short c;
            if (laplace.locate(events[i], coord))
            {
                c = laplace.classify(coord);
                if (c > results->valid_clusters)
                    c = 0;
            }
            else
                c = results->valid_clusters + 1;
            results->cluster_events[c].push_back(i);
        }
        for (unsigned c = 0; c < results->cluster_events.size(); ++c)
        {
            if (results->cluster_events.size() < selections.min_events)
            {
                for (unsigned i : results->cluster_events[c])
                    results->cluster_events[0].push_back(i);
                results->cluster_events[c].clear();
                results->valid_clusters--;
            }
        }

        for_each_plane<Dimension>([this, results, data, pop_name](unsigned i, unsigned j)
        {
            results->future_marginals.push_back(compute_plane.enqueue([this, results, data, pop_name, i, j]() { return do_Marginal(results, data, i, j, pop_name); }));
        });

        for (unsigned c = 0; c <= results->valid_clusters; ++c)
            if (results->cluster_events[c].size() > 0)
                for (unsigned i = 0; i < Dimension; i++)
                    for (unsigned e : results->cluster_events[c])
                        results->means[c][i] += events[e][i];
        for (unsigned c = 0; c <= results->valid_clusters; ++c)
            if (results->cluster_events[c].size() > 0)
                for (unsigned i = 0; i < Dimension; i++)
                    results->means[c][i] /= results->cluster_events[c].size();
        for (unsigned c = 0; c <= results->valid_clusters; ++c)
            if (results->cluster_events[c].size() > 1)
                for (unsigned i = 0; i < Dimension; i++)
                    for (unsigned j = 0; j < Dimension; j++)
                        for (unsigned e : results->cluster_events[c])
                            results->covariances[c][i][j] += (events[e][i] - results->means[c][i]) * (events[e][j] - results->means[c][j]);
        for (unsigned c = 0; c <= results->valid_clusters; ++c)
            if (results->cluster_events[c].size() > 1)
                for (unsigned i = 0; i < Dimension; i++)
                    for (unsigned j = 0; j < Dimension; j++)
                        results->covariances[c][i][j] /= results->cluster_events[c].size() - 1;


        std::cout << "Found " << results->valid_clusters << " valid clusters." << std::endl;
        return results;
    }
    
    int run();
};
