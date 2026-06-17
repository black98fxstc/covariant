#pragma once

#include <vector>
#include <string>
#include <memory>
#include <future>
#include <atomic>
#include <limits>
#include <iostream>
#include <cstdint>

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

    void wait() noexcept;
};

template <unsigned Dimension>
void for_each_plane(std::function<void(const unsigned i, const unsigned j)> func)
{
    for (unsigned i = 0; i < Dimension - 1; i++)
        for (unsigned j = i + 1; j < Dimension; j++)
            func(i, j);
};

class Marginal_Results
{
public:
    std::vector<std::vector<std::vector<double>>> class_data;
    std::vector<std::vector<double>> quant_data;
    std::string x_label, y_label;
    unsigned clusters;
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
    ThreadPool compute_plane{1};//std::thread::hardware_concurrency()};
    ThreadPool control_plane{4};

    int parse_args(int argc, char *argv[]);

    Qualify_Results do_Qualify(const std::vector<float> *data, const Measurement X, const std::vector<bool> &included, std::string pop_name);

    Projection_Results do_Projection(const std::vector<std::vector<float> *> &data, const Measurement X, const Measurement Y, const std::vector<bool> &included, std::string pop_name);

    Marginal_Results do_Marginal(const std::vector<std::vector<float> *> &data, const std::string pop_name) noexcept
    {
        Marginal_Results results;
        return results;
    }

    Pursuit_Results do_Pursuit(const std::vector<std::vector<float> *> &data, const std::vector<bool> &included, std::string pop_name);

    template <unsigned Dimension>
    Laplace_Results do_Laplace(const std::vector<std::vector<float> *> &data, const std::vector<bool> &included, std::string pop_name)
    {
        Laplace_Results results;
        results.parent_name = pop_name;
        results.idx.reserve(included.size());
        for (size_t i = 0; i < included.size(); ++i)
            if (included[i])
                results.idx.push_back(i);

        if (results.idx.empty()) {
            std::cout << "No valid events found." << std::endl;
            return results;
        }

        std::cout << "Begin Laplacian clustering on " << results.idx.size() << " events..." << std::endl;
        Events<Dimension> events;
        events.resize(results.idx.size());
        for (unsigned d = 0; d < Dimension; ++d)
            for (size_t i = 0; i < events.size(); ++i)
                events[i][d] = (*data[d])[results.idx[i]];

        Laplace<Dimension> laplace(selections.grid_size);
        laplace.verify = params.verify;
        laplace.antialias = params.antialias;

        size_t valid = 0;
        for (const auto &e : events)
            if (laplace.event(e))
                valid++;
        std::cout << "Found " << valid << " valid events..." << std::endl;

        std::cout << "Calculating the Laplacian of the sample..." << std::endl;
        laplace.analyze(selections.smoothing, selections.threshold);
        if (laplace.differentialError() > .0001)
            std::cout << "Differential equation solution is unusually bad " << laplace.differentialError() << std::endl;
        else
            std::cout << "Consistency checks passed..." << std::endl;

        std::cout << "Performing Laplacian clustering..." << std::endl;
        results.clusters_found = laplace.cluster(selections.threshold);
        results.valid_clusters = std::min(selections.max_clusters, results.clusters_found);
        results.cluster_events.resize(results.valid_clusters + 2);
        results.means.resize(results.valid_clusters + 1);
        results.covariances.resize(results.valid_clusters + 1);
        for (unsigned c = 0; c <= results.valid_clusters; ++c)
        {
            results.means[c].resize(Dimension);
            results.covariances[c].resize(Dimension);
            for (unsigned i = 0; i < Dimension; i++)
                    results.covariances[c][i].resize(Dimension);
        }
        
        Coordinates coord(laplace);
        for (size_t i = 0; i < results.idx.size(); ++i)
        {
            unsigned short c;
            if (laplace.locate(events[i], coord))
            {
                c = laplace.classify(coord);
                if (c > results.valid_clusters)
                    c = 0;
            }
            else
                c = results.valid_clusters + 1;
            results.cluster_events[c].push_back(i);
        }
        for (unsigned c = 0; c < results.cluster_events.size(); ++c)
        {
            if (results.cluster_events.size() < selections.min_events)
            {
                for (unsigned i : results.cluster_events[c])
                    results.cluster_events[0].push_back(i);
                results.cluster_events[c].clear();
                results.valid_clusters--;
            }
        }

        for_each_plane<Dimension>([this, &results, &data, &pop_name](unsigned i, unsigned j)
        {
            results.future_marginals.push_back(compute_plane.enqueue([this, data, pop_name]() { return do_Marginal(data, pop_name); }));
        });
        
        for (unsigned c = 0; c <= results.valid_clusters; ++c)
            if (results.cluster_events[c].size() > 0)
                for (unsigned i = 0; i < Dimension; i++)
                    for (unsigned e : results.cluster_events[c])
                        results.means[c][i] += events[e][i];
        for (unsigned c = 0; c <= results.valid_clusters; ++c)
            if (results.cluster_events[c].size() > 0)
                for (unsigned i = 0; i < Dimension; i++)
                    results.means[c][i] /= results.cluster_events[c].size();
        for (unsigned c = 0; c <= results.valid_clusters; ++c)
            if (results.cluster_events[c].size() > 1)
                for (unsigned i = 0; i < Dimension; i++)
                    for (unsigned j = 0; j < Dimension; j++)
                        for (unsigned e : results.cluster_events[c])
                            results.covariances[c][i][j] += (events[e][i] - results.means[c][i]) * (events[e][j] - results.means[c][j]);
        for (unsigned c = 0; c <= results.valid_clusters; ++c)
            if (results.cluster_events[c].size() > 1)
                for (unsigned i = 0; i < Dimension; i++)
                    for (unsigned j = 0; j < Dimension; j++)
                        results.covariances[c][i][j] /= results.cluster_events.size() - 1;


        std::cout << "Found " << results.valid_clusters << " valid clusters." << std::endl;
        return results;
    }
    
    int run();
};
