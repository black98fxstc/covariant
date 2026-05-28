#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <cxxopts.hpp>
#include <filesystem>
#include <iomanip>
#include <nlohmann/json.hpp>
#include <cstdlib>
#include <matplot/matplot.h>

#include "Covariant.hpp"

using json = nlohmann::json;

namespace fs = std::filesystem;

struct Params
{
    std::string filename;
    std::string out_dir;
    unsigned dimension;
    unsigned grid;
    unsigned points;
    float smooth;
    float threshold;
    bool visual;
    unsigned max_clusters;
    size_t min_events;
    bool verbose;
    bool antialias;
    bool verify;
    bool grow;
    bool ascii;
};

template <unsigned Dimension>
void for_each_plane(std::function<void(const unsigned i, const unsigned j)> func)
{
    for (unsigned i = 0; i < Dimension - 1; i++)
        for (unsigned j = i + 1; j < Dimension; j++)
            func(i, j);
};

// HSL (Hue [0,1], Saturation [0,1], Lightness [0,1]) to RGB [0,1]
std::array<double, 3> hsl2rgb(double h, double s, double l)
{
    double c = (1.0 - std::abs(2.0 * l - 1.0)) * s;
    double hp = h * 6.0;
    double x = c * (1.0 - std::abs(std::fmod(hp, 2.0) - 1.0));
    double r1 = 0, g1 = 0, b1 = 0;
    if (hp < 1)
    {
        r1 = c;
        g1 = x;
    }
    else if (hp < 2)
    {
        r1 = x;
        g1 = c;
    }
    else if (hp < 3)
    {
        g1 = c;
        b1 = x;
    }
    else if (hp < 4)
    {
        g1 = x;
        b1 = c;
    }
    else if (hp < 5)
    {
        r1 = x;
        b1 = c;
    }
    else
    {
        r1 = c;
        b1 = x;
    }
    double m = l - c / 2.0;
    return {r1 + m, g1 + m, b1 + m};
}

// RGB [0,1] to HSL (Hue [0,1], Saturation [0,1], Lightness [0,1])
std::array<double, 3> rgb2hsl(double r, double g, double b)
{
    double max = std::max({r, g, b});
    double min = std::min({r, g, b});
    double h, s, l = (max + min) / 2.0;
    if (max == min)
        h = s = 0;
    else
    {
        double d = max - min;
        s = l > 0.5 ? d / (2.0 - max - min) : d / (max + min);
        if (max == r)
            h = (g - b) / d + (g < b ? 6 : 0);
        else if (max == g)
            h = (b - r) / d + 2;
        else
            h = (r - g) / d + 4;
        h /= 6.0;
    }
    return {h, s, l};
}

void make_plot(const std::string &path, const std::vector<std::vector<double>> &class_data, const std::vector<std::vector<double>> &quant_data)
{
    using namespace matplot;

    auto fig = figure(false); // Off-screen rendering
    fig->size(800, 800);
    auto ax = fig->current_axes();

    // Map the classification matrix to the [0, 1] range for both axes.
    // Using the global function with the axes handle resolves AxesHandle ambiguities.
    imagesc(ax, 0.0, 1.0, 0.0, 1.0, class_data);
    ax->hold(on);

    // Generate coordinate matrices for contour mapping.
    auto x_range = linspace(0, 1, class_data[0].size());
    auto y_range = linspace(0, 1, class_data.size());
    auto [X, Y] = meshgrid(x_range, y_range);

    // Overlay contours.
    auto c = ax->contour(X, Y, quant_data);
    c->color("black");
    c->line_width(1.0);

    // Set explicit limits for the axes
    ax->xlim({0, 1});
    ax->ylim({0, 1});

    fig->save(path);
}

template <unsigned Dimension>
int do_it(const Params &params)
{
    typename Weighty<Dimension>::Events events;
    std::string ext = params.ascii ? ".txt" : ".dat";
    std::cout << "Loading events from " << params.filename << ext << "..." << std::endl;
    if (!events.read(params.filename + ext, params.ascii))
    {
        std::cerr << "Error: Could not open event file: " << params.filename << ext << std::endl;
        return 1;
    }
    std::cout << "Loaded " << events.size() << " events." << std::endl;

    // Create output directory
    fs::create_directories(params.out_dir);

    // 1. Initial Riemann object to perform global clustering via Laplace
    Riemann<Dimension> global(params.grid);
    global.visualize = false; // We will save files manually into the dir
    global.antialias = params.antialias;
    global.verify = params.verify;

    for (const auto &e : events)
    {
        global.event(e);
    }

    std::cout << "Initial global analysis for clustering..." << std::endl;
    global.Laplace<Dimension>::analyze(params.smooth, params.threshold);
    unsigned num_clusters = global.cluster(params.threshold, params.grow);
    std::cout << "Found " << num_clusters << " clusters. Starting per-cluster analysis..." << std::endl;

    // Save the global classification volume once
    global.cluster_id.write(params.out_dir + "/" + params.filename + "_classes.bin");

    json report;
    report["filename"] = params.filename;
    report["dimension"] = Dimension;
    report["grid_size"] = params.grid;
    report["total_events"] = events.size();
    report["num_clusters"] = num_clusters;
    report["clusters"] = json::array();

    std::ofstream xml_out(params.filename + ".report.xml");
    xml_out << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
    xml_out << "<DarwinReport file=\"" << params.filename << "\" dimensions=\"" << Dimension << "\">\n";
    xml_out << "  <Summary totalEvents=\"" << events.size() << "\" clustersFound=\"" << num_clusters << "\" />\n";
    xml_out << "  <Clusters>\n";

    // Prepare HTML report
    std::ofstream html_out(params.out_dir + "/index.html");
    html_out << "<html><head><title>Darwin Report: " << params.filename << "</title>";
    html_out << "<style>body{font-family:sans-serif;margin:40px;} table{border-collapse:collapse;width:100%;} th,td{padding:10px;border:1px solid #ccc;text-align:left;} img{max-width:800px;}</style></head>";
    html_out << "<body><h1>Darwin Report: " << params.filename << "</h1>";
    html_out << "<table><tr><th>ID</th><th>Events</th><th>Diff Error</th><th>Factor Prob</th><th>Mean Vector</th><th>Visualization</th></tr>";

    // Setup for 2D projections
    Weighty<2> marginal(params.grid);
    // typename Weighty<2>::Event marginal_event;
    marginal.visualize = false; // We will save files manually into the dir
    marginal.antialias = params.antialias;
    marginal.verify = params.verify;

    std::string proj_dir = params.out_dir + "/projections";
    fs::create_directories(proj_dir);
    // unsigned grid_pts = grid_size + 1;

    std::vector<std::vector<typename Weighty<Dimension>::Event>> cluster_events(num_clusters + 1);
    for (const auto &e : events)
    {
        unsigned short c = 0;
        Coordinate<Dimension> coord(global);
        if (global.locate(e, coord))
        {
            size_t idx = (size_t)coord;
            c = static_cast<const Function<Dimension, unsigned short> &>(global.cluster_id)[idx];
        }
        else
            c = 0;

        if (c < num_clusters + 1)
            cluster_events[c].push_back(e);
        else
            cluster_events[0].push_back(e);
    }
    for_each_plane<Dimension>([&params, &num_clusters, &cluster_events, &marginal, &proj_dir](unsigned i, unsigned j)
                              {
        marginal.reset();
        Coordinate<2> marginal_coord(marginal);
        std::vector<unsigned short> marginal_klass(params.points * params.points, 0);

        for (unsigned short c = 0; c <= num_clusters; ++c)
        {
            Weighty<2>::Event marginal_event;
            for (auto &e : cluster_events[c])
            {
                marginal_event[0] = e[i];
                marginal_event[1] = e[j];
                marginal.event(marginal_event);
                if (c == 0) continue;

                if (marginal.locate(marginal_event, marginal_coord))
                {
                    size_t idx = (size_t)marginal_coord;
                    unsigned short d = marginal_klass[idx];
                    if (d == 0) marginal_klass[idx] = c;
                    else marginal_klass[idx] = std::min(c, d);
                }
            }
        }
        marginal.prepare(params.smooth);

        if (params.visual) {
            std::vector<std::vector<double>> class_data(params.points, std::vector<double>(params.points));
            std::vector<std::vector<double>> quant_data(params.points, std::vector<double>(params.points));
            for (unsigned y = 0; y < params.points; ++y) {
                for (unsigned x = 0; x < params.points; ++x) {
                    size_t idx = x + y * params.points;
                    class_data[y][x] = (double)marginal_klass[idx];
                    quant_data[y][x] = (double)static_cast<const Function<2, float>&>(marginal.quantile)[idx];
                }
            }

            std::string path = proj_dir + "/full_d" + std::to_string(i) + "_d" + std::to_string(j) + ".png";
            make_plot(path, class_data, quant_data);
        } });
    for (unsigned c = 1; c <= std::min(num_clusters, params.max_clusters); ++c)
    {
        size_t cluster_event_count = cluster_events[c].size();
        if (cluster_events[c].size() < params.min_events)
        {
            std::cout << "  Skipping small " << cluster_event_count << "events  cluster " << c << "..." << std::endl;
            continue;
        }

        // Riemannian analysis on the cluster subset
        global.reset();

        std::vector<double> mean(Dimension, 0.0);
        std::vector<std::vector<double>> covariance(Dimension, std::vector<double>(Dimension, 0.0));
        for (auto &e : cluster_events[c])
        {
            global.event(e);
            for (unsigned i = 0; i < Dimension; i++)
                mean[i] += e[i];
        }
        for (unsigned i = 0; i < Dimension; i++)
            mean[i] /= cluster_events[c].size();
        for (auto &e : cluster_events[c])
            for (unsigned i = 0; i < Dimension; i++)
                for (unsigned j = 0; j < Dimension; j++)
                    covariance[i][j] += (e[i] - mean[i]) * (e[j] - mean[j]);
        for (unsigned i = 0; i < Dimension; i++)
            for (unsigned j = 0; j < Dimension; j++)
                covariance[i][j] /= (cluster_events[c].size() - 1);

        global.Riemann<Dimension>::analyze(params.smooth, params.threshold);

        for_each_plane<Dimension>([&params, c, &cluster_events, &marginal, &proj_dir](unsigned i, unsigned j)
                                  {
            Weighty<2>::Event marginal_event;
            Coordinate<2> marginal_coord(marginal);
            marginal.reset();
            for (auto &e : cluster_events[c])
            {
                marginal_event[0] = e[i];
                marginal_event[1] = e[j];
                marginal.event(marginal_event);
            }
            marginal.prepare(params.smooth);

            // make cluster marginal plot
            if (params.visual) {
                std::vector<std::vector<double>> class_data(params.points, std::vector<double>(params.points));
                std::vector<std::vector<double>> quant_data(params.points, std::vector<double>(params.points));

                for (unsigned y = 0; y < params.points; ++y) {
                    for (unsigned x = 0; x < params.points; ++x) {
                        size_t idx = x + y * params.points;
                        class_data[y][x] = (double)c;
                        quant_data[y][x] = (double)static_cast<const Function<2, float>&>(marginal.quantile)[idx];
                    }
                }

                std::string path = proj_dir + "/cluster" + std::to_string(c) + "_d" + std::to_string(i) + "_d" + std::to_string(j) + ".png";
                make_plot(path, class_data, quant_data);
            } });

        std::string img_rel_name = params.filename + ".cluster" + std::to_string(c) + ".png";
        std::string img_full_path = params.out_dir + "/" + img_rel_name;

        // Populate JSON report
        json c_info;
        c_info["id"] = c;
        c_info["event_count"] = cluster_event_count;
        c_info["differential_error"] = global.differentialError();
        c_info["factor_probability"] = global.factorProbability();
        c_info["mean"] = mean;
        c_info["covariance"] = covariance;
        c_info["image"] = img_rel_name;
        report["clusters"].push_back(c_info);

        // Populate XML data
        xml_out << "    <Cluster id=\"" << c << "\" events=\"" << cluster_event_count << "\">\n";
        xml_out << "      <Metrics diffError=\"" << global.differentialError() << "\" ";
        xml_out << "factorProb=\"" << global.factorProbability() << "\" />\n";
        xml_out << "      <Mean>";
        for (unsigned i = 0; i < Dimension; ++i)
            xml_out << (i == 0 ? "" : " ") << mean[i];
        xml_out << "</Mean>\n";
        xml_out << "      <Visualization src=\"" << img_rel_name << "\" />\n";
        xml_out << "    </Cluster>\n";

        // Update HTML
        html_out << "<tr><td>" << c << "</td><td>" << cluster_event_count << "</td>";
        html_out << "<td>" << global.differentialError() << "</td><td>" << global.factorProbability() << "</td>";
        html_out << "<td>";
        for (unsigned i = 0; i < Dimension; ++i)
        {
            html_out << (i == 0 ? "" : ", ") << std::fixed << std::setprecision(3) << mean[i];
        }
        html_out << "</td>";
        html_out << "<td><img src=\"" << img_rel_name << "\"></td></tr>";

        if (params.visual)
        {
            std::cout << "      Generating visualization: " << img_rel_name << "..." << std::endl;
            // Pass parameters to MATLAB script
            std::string cmd = "matlab -batch \"plot_covariant_projections('" + params.out_dir + "', '" + params.filename + "', " +
                              std::to_string(c) + ", " + std::to_string(Dimension) + ", '" + img_rel_name + "');\"";
            std::system(cmd.c_str());
        }
    }

    html_out << "</table></body></html>";
    html_out.close();

    xml_out << "  </Clusters>\n";
    xml_out << "</DarwinReport>\n";
    xml_out.close();

    std::ofstream json_out(params.filename + ".report.json");
    json_out << report.dump(4);
    json_out.close();

    std::cout << "Analysis complete. Reports saved." << std::endl;

    return 0;
}

int main(int argc, char *argv[])
{
    Params params;

    cxxopts::Options options("Darwin", "Hierarchical Laplacian and Riemannian analysis");

    options.add_options()("f,file", "Input filename", cxxopts::value<std::string>()->default_value("test_data"))("d,dimension", "Dimension", cxxopts::value<unsigned>()->default_value("2"))("g,grid", "Grid resolution", cxxopts::value<unsigned>())("s,smooth", "Smoothing factor", cxxopts::value<float>()->default_value("0.01"))("t,threshold", "Threshold", cxxopts::value<float>()->default_value("0.001"))("visual", "Save visualization files", cxxopts::value<bool>()->default_value("false"))("max-clusters", "Maximum number of clusters to report", cxxopts::value<unsigned>()->default_value("10"))("min-events", "Minimum number of events for a cluster to be reported", cxxopts::value<size_t>()->default_value("1"))("v,verbose", "Verbose", cxxopts::value<bool>()->default_value("false"))("antialias", "Antialiasing", cxxopts::value<bool>()->default_value("true"))("verify", "Verify consistency", cxxopts::value<bool>()->default_value("true"))("grow", "Grow clusters", cxxopts::value<bool>()->default_value("false"))("a,ascii", "Use ASCII data", cxxopts::value<bool>()->default_value("false"))("h,help", "Print usage");

    options.parse_positional({"file"});

    auto result = options.parse(argc, argv);

    if (result.count("help"))
    {
        std::cout << options.help() << std::endl;
        return 0;
    }

    params.filename = result["file"].as<std::string>();
    params.dimension = result["dimension"].as<unsigned>();
    params.smooth = result["smooth"].as<float>();
    params.threshold = result["threshold"].as<float>();
    params.visual = result["visual"].as<bool>();
    params.max_clusters = result["max-clusters"].as<unsigned>();
    params.min_events = result["min-events"].as<size_t>();
    params.antialias = result["antialias"].as<bool>();
    params.verify = result["verify"].as<bool>();
    params.grow = result["grow"].as<bool>();
    params.verbose = result["verbose"].as<bool>();
    params.ascii = result["ascii"].as<bool>();

    params.out_dir = params.filename + ".darwin";
    if (result.count("grid"))
        params.grid = result["grid"].as<unsigned>();
    else
        switch (params.dimension)
        {
        case 2:
            params.grid = 256;
            break;
        case 3:
            params.grid = 64;
            break;
        default:
            params.grid = 32;
            break;
        }
    params.points = params.grid + 1;

    std::cout << "Darwin running: " << params.filename << " (Dim: " << params.dimension << ")" << std::endl;
    std::cout << "  Max clusters to report: " << params.max_clusters << std::endl;
    std::cout << "  Min events per cluster: " << params.min_events << std::endl;

    switch (params.dimension)
    {
    case 2:
        return do_it<2>(params);
    case 3:
        return do_it<3>(params);
    case 4:
        return do_it<4>(params);
    default:
        std::cerr << "Unsupported dimension." << std::endl;
        return 1;
    }
}
// ```

// ### Summary of Changes:
// 1.  **Makefile**: Added the `darwin` target and included it in the `all` and `clean` rules.
// 2.  **darwin.cpp**:
//    *   It starts by loading the event data.
//    *   It creates a global `Riemann` object to perform the initial `analyze()` and `cluster()` calls.
//    *   It iterates through each discovered cluster, filters the original events to find those belonging to that cluster, and runs a secondary `Riemann::analyze()` on that subset.
//    *   Results are collected into a `nlohmann::json` object and written to a `.json` file.
//    *   A simple XML report is also generated, which can be easily rendered by a browser using a basic CSS or XSLT stylesheet.

// <!--
// [PROMPT_SUGGESTION]Can you add an XSLT stylesheet to the XML output so it renders as a table in a web browser?[/PROMPT_SUGGESTION]
// [PROMPT_SUGGESTION]How can I optimize the cluster filtering loop in darwin.cpp to avoid iterating over all events for every cluster?[/PROMPT_SUGGESTION]
