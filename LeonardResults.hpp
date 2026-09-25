#pragma once

#include <cstddef>
#include <fstream>
#include <future>
#include <limits>
#include <memory>
#include <mutex>
#include <string>
#include <utility>
#include <vector>

#include <nlohmann/json.hpp>

#include "Geometry.hpp"

using json = nlohmann::json;

class Qualify_Results
{
public:
    Measurement X;
    double KLDn = 0;
    double KLDe = 0;
    bool qualified = false;

    Qualify_Results() = default;
    Qualify_Results(const Measurement X) noexcept : X(X) {}
};

class EPP_Node_Results
{
public:
    std::string image_in;
    std::string image_out;
    std::vector<std::future<void>> future_plots;

    void wait_for_plots();

    EPP_Node_Results() = default;
    EPP_Node_Results(EPP_Node_Results &&) = default;
    EPP_Node_Results &operator=(EPP_Node_Results &&) = default;
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
    Measurement X, Y;
    Polygon separatrix;
    double score = std::numeric_limits<double>::infinity();
    struct Gating
    {
        size_t count = 0;
        std::vector<bool> set;
        Polygon polygon;
    } in, out;
    unsigned int pass = 0, clusters = 0, graphs = 0, merges = 0, splits = 0;

    Projection_Results() = default;
    Projection_Results(Measurement X, Measurement Y) noexcept
        : X(X < Y ? X : Y), Y(X < Y ? Y : X) {}

private:
    void close_clockwise(Polygon &polygon) const noexcept;
};

class Pursuit_Results
{
public:
    std::string node_id = "1";
    std::string branch = "root";
    size_t event_count = 0;
    size_t total_events = 0;
    std::vector<double> means;
    bool is_leaf = true;
    Measurement gate_x = 0;
    Measurement gate_y = 0;
    bool has_gate = false;
    std::string polygon_image;
    double pct_parent = 0, pct_total = 0;

    std::vector<Measurement> qualified;
    double best_score = 0;
    std::unique_ptr<Projection_Results> best_split;
    std::vector<std::future<Pursuit_Results>> future_children;
    std::vector<Pursuit_Results> children;
    std::future<EPP_Node_Results> future_node;
    std::unique_ptr<EPP_Node_Results> EPP_node;
    std::vector<std::future<void>> future_plots;
    std::vector<std::string> sample_images;

    void wait_for_results() noexcept;
    void wait_for_plots();
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

    Marginal_Results() = default;
    Marginal_Results(std::shared_ptr<Laplace_Results> laplace) noexcept : laplace(laplace) {}
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
    std::vector<std::future<void>> future_plots;
    std::vector<std::string> sample_images;
    std::vector<std::string> cluster_images;
    std::shared_ptr<std::mutex> image_mutex = std::make_shared<std::mutex>();

    void add_sample_image(std::string image)
    {
        std::lock_guard<std::mutex> lock(*image_mutex);
        sample_images.push_back(std::move(image));
    }

    void add_cluster_image(std::string image)
    {
        std::lock_guard<std::mutex> lock(*image_mutex);
        cluster_images.push_back(std::move(image));
    }

    void wait_for_results() noexcept
    {
        for (auto &future : future_marginals)
            marginals.push_back(future.get());
    }

    void wait_for_plots()
    {
        for (auto &future : future_plots)
            if (future.valid())
                future.get();
    }

    Laplace_Results() = default;
    Laplace_Results(Laplace_Results &&) = default;
    Laplace_Results &operator=(Laplace_Results &&) = default;
};
