#pragma once

#include <algorithm>
#include <atomic>
#include <cstddef>
#include <functional>
#include <memory>
#include <string>
#include <thread>
#include <vector>

#include "Geometry.hpp"
#include "LeonardResults.hpp"
#include "FlowJo.hpp"
#include "Workers.hpp"
#include "Covariant.hpp"
#include "Events.hpp"
#include "Log.hpp"

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
    std::string img_dir = "images";
};

// Forward declare plotting functions implemented in Reports.cpp
void make_marginal_plot(const std::string &path, const std::vector<std::vector<std::vector<double>>> &class_data, const std::vector<std::vector<double>> &quant_data);
void make_overlay_transparent(const std::string &path);
void write_rgb_png(const std::string &path, const std::vector<std::vector<std::vector<double>>> &class_data);
void make_gating_plot(const std::string &path, const std::vector<std::vector<double>> &quant_data, const Measurement X, const Measurement Y, const Polygon &polygon);

template <unsigned Dimension>
void for_each_plane(std::function<void(const unsigned i, const unsigned j)> func)
{
    for (unsigned i = 0; i < Dimension - 1; i++)
        for (unsigned j = i + 1; j < Dimension; j++)
            func(i, j);
};

class Leonard
{
    // Declared first so it outlives the thread pools, which can still log while they shut down.
    SessionLog session_;

public:
    Leonard() : say(&session_, LogLevel::Say), log(&session_, LogLevel::Log) {}
    Leonard(const Leonard &) = delete;
    Leonard &operator=(const Leonard &) = delete;

    // say goes to the session file and to the console unless --quiet.
    // log goes to the session file and to the console only with --verbose.
    LogChannel say;
    LogChannel log;

    Params params;
    SelectionState selections;
    size_t laplacian_offset = 0;
    std::atomic<unsigned> laplace_id_counter = 0;
    std::atomic<unsigned> EPP_id_counter = 0;
    Workspace ws;
    std::vector<SampleData> dummy_samples;
    std::vector<std::vector<double>> colors;
    ThreadPool compute_plane{std::thread::hardware_concurrency()};
    ThreadPool control_plane{4};
    ThreadPool plot_plane{1}; // std::max(1u, std::thread::hardware_concurrency())};

    int parse_args(int argc, char *argv[]);

    Qualify_Results do_Qualify(const std::vector<float> *data, const Measurement X, const std::vector<bool> &included, std::string pop_name);

    Projection_Results do_Projection(const std::vector<std::vector<float> *> &data, const Measurement X, const Measurement Y, const std::vector<bool> &included, std::string pop_name);

    EPP_Node_Results do_EPP_Node(const std::vector<std::vector<float> *> data, const std::vector<bool> &included, const Measurement X, const Measurement Y, const Polygon &in_poly, const Polygon &out_poly, std::string pop_name, std::string node_id);

    Pursuit_Results do_Pursuit(const std::vector<std::vector<float> *> &data, std::vector<bool> included, std::string pop_name, size_t total_events, std::string node_id = "1", std::string branch = "root");

    Marginal_Results do_Marginal(std::shared_ptr<Laplace_Results> laplace, const std::vector<std::vector<float> *> &data, const Measurement i, const Measurement j, std::string pop_name) noexcept;

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

        say << "Begin Laplacian clustering on " << results->idx.size() << " events..." << std::endl;
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
        say << "Found " << valid_events << " valid events..." << std::endl;

        say << "Calculating the Laplacian of the sample..." << std::endl;
        laplace.analyze(selections.smoothing, selections.threshold);
        if (laplace.differentialError() > .0001)
            say << "Differential equation solution is unusually bad " << laplace.differentialError() << std::endl;
        else
            say << "Consistency checks passed..." << std::endl;

        say << "Performing Laplacian clustering..." << std::endl;
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
            if (results->cluster_events[c].size() < selections.min_events)
            {
                for (unsigned i : results->cluster_events[c])
                    results->cluster_events[0].push_back(i);
                results->cluster_events[c].clear();
                results->valid_clusters--;
            }
        }
        results->classification.assign(results->idx.size(), 0);
        for (unsigned c = 0; c < results->cluster_events.size(); ++c)
            for (unsigned local : results->cluster_events[c])
                results->classification[local] = static_cast<unsigned short>(c);

        for_each_plane<Dimension>([this, results, data, pop_name](unsigned i, unsigned j)
                                  { results->future_marginals.push_back(compute_plane.enqueue([this, results, data, pop_name, i, j]()
                                                                                              { return do_Marginal(results, data, i, j, pop_name); })); });

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

        say << "Found " << results->valid_clusters << " valid clusters." << std::endl;
        return results;
    }

    int run();
};
