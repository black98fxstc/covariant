#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <cxxopts.hpp>
#include <filesystem>
#include <iomanip>
#include <nlohmann/json.hpp>
#include <cstdlib>
#include <unistd.h>
#include <fcntl.h>
#include <matplot/matplot.h>
#include <sys/wait.h>
#include <thread>
#include <algorithm>

#include "Covariant.hpp"

using json = nlohmann::json;

namespace fs = std::filesystem;

struct Params
{
    std::string filename;
    std::string out_dir;
    std::string log_path;
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
    unsigned num_clusters = 0;
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

/**
 * @brief Simple HSV to RGB conversion helper.
 * 
 * @param h Hue [0, 1]
 * @param s Saturation [0, 1]
 * @param v Value [0, 1]
 * @return std::vector<double> RGB components
 */
std::vector<double> hsv_to_rgb(double h, double s, double v) {
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

/**
 * @brief Generates a cluster-aware colormap similar to MATLAB's [1 1 1; hueMap(n)].
 * 
 * @param num_clusters Number of actual clusters found in the data.
 * @param max_n The maximum capacity for the hue map (e.g., 10 for qmap).
 * @return std::vector<std::vector<double>> The colormap matrix.
 */
std::vector<std::vector<double>> get_cluster_qmap(const Params &params) {
    // Determine how many hues we actually need to generate
    size_t n = std::min((size_t)params.num_clusters, (size_t)params.max_clusters);
    
    std::vector<std::vector<double>> cmap;
    cmap.reserve(n + 1);

    // 1. First element is always white (index 0 / Background / Unassigned)
    cmap.push_back({1.0, 1.0, 1.0});

    // 2. Generate the shifted hues (equivalent to hueMap(n))
    // We replicate linspace(0.6666, 0, n)
    for (size_t i = 0; i < n; ++i) {
        double h = 0.6666;
        if (n > 1) {
            h = 0.6666 * (1.0 - static_cast<double>(i) / (n - 1));
        }
        cmap.push_back(hsv_to_rgb(h, 1.0, 0.75));
    }

    return cmap;
}

// RAII helper to silence Gnuplot terminal noise while capturing logs/errors in a file.
struct GnuplotSilencer
{
    int old_stdout = -1;
    int old_stderr = -1;
    int dev_null = -1;
    int log_fd = -1;

    GnuplotSilencer(const std::string &log_path)
    {
        std::cout.flush();
        std::cerr.flush();
        fflush(stdout);
        fflush(stderr);

        old_stdout = dup(STDOUT_FILENO);
        old_stderr = dup(STDERR_FILENO);
        dev_null = open("/dev/null", O_WRONLY);
        log_fd = open(log_path.c_str(), O_WRONLY | O_CREAT | O_APPEND, 0644);

        if (dev_null != -1 && old_stdout != -1) dup2(dev_null, STDOUT_FILENO);
        if (log_fd != -1 && old_stderr != -1) dup2(log_fd, STDERR_FILENO);
    }

    ~GnuplotSilencer()
    {
        std::cout.flush();
        std::cerr.flush();
        fflush(stdout);
        fflush(stderr);

        if (old_stdout != -1) { dup2(old_stdout, STDOUT_FILENO); close(old_stdout); }
        if (old_stderr != -1) { dup2(old_stderr, STDERR_FILENO); close(old_stderr); }
        if (dev_null != -1) close(dev_null);
        if (log_fd != -1) close(log_fd);
    }
};

void make_plot(const std::string &path, const Params &params, const std::vector<std::vector<double>> &class_data, const std::vector<std::vector<double>> &quant_data)
{
    using namespace matplot;

    // Instantiate the silencer to capture Gnuplot stderr in the log file
    // GnuplotSilencer silencer(params.log_path);

    auto fig = figure(true);   // Quiet mode    // fig->visible(false);       // Explicitly disable windowing for headless/CLI environments
    fig->size(400, 400);
    auto ax = fig->current_axes();
    auto qmap = get_cluster_qmap(params);
    ax->colormap(qmap);
    ax->color_box_range(0.0, static_cast<double>(qmap.size()));
    ax->max_colors(qmap.size() - 1);

    // Enable hold before plotting to buffer commands and prevent intermediate 
    // output from leaking into the terminal during figure initialization.
    hold(ax, on);

    // Map the classification matrix to the [0, 1] range for both axes.
    // Using the global function with the axes handle resolves AxesHandle ambiguities.
    imagesc(ax, 0.0, 1.0, 0.0, 1.0, class_data);
    // Ensure cluster IDs map to the correct colormap indices

    // Generate coordinate matrices for contour mapping.
    auto x_range = linspace(0, 1, class_data[0].size());
    auto y_range = linspace(0, 1, class_data.size());
    auto [X, Y] = meshgrid(x_range, y_range);

    // Overlay contours.
    auto c = contour(ax, X, Y, quant_data);
    c->levels(matplot::iota(0.1, 0.1, 0.9)); // Equivalent to 0.1:0.1:0.9
    c->color("black");
    c->line_width(1.2);

    // auto c2 = contour(ax, X, Y, quant_data);
    // c2->levels({0.001});
    // c2->color("red");
    // c2->line_width(1.2);

    // auto c3 = contour(ax, X, Y, quant_data);
    // c3->levels({0.01});
    // c3->color("red");
    // c3->line_width(1.2);

    // Set explicit limits for the axes
    ax->xlim({0, 1});
    ax->ylim({0, 1});

    fig->save(path);
}

template <unsigned Dimension>
int do_it(Params &params)
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

    // Create logs folder and truncate the session log file
    std::string log_dir = params.out_dir + "/logs";
    fs::create_directories(log_dir);
    params.log_path = log_dir + "/gnuplot.log";
    std::ofstream(params.log_path, std::ios::trunc);

    // Riemann object for the whole n-dimensional sample
    Riemann<Dimension> global(params.grid);
    global.visualize = false; // We will save files manually into the dir
    global.antialias = params.antialias;
    global.verify = params.verify;

    // Weighty object for 2-dimensional marginal distributions
    Weighty<2> marginal(params.grid);
    marginal.visualize = false; // We will save files manually into the dir
    marginal.antialias = params.antialias;
    marginal.verify = params.verify;

    size_t valid_events = 0;
    for (const auto &e : events)
        if (global.event(e)) ++valid_events;

    std::cout << "Initial global analysis found " << valid_events << " valid events..." << std::endl;
    global.Laplace<Dimension>::analyze(params.smooth, params.threshold);
    params.num_clusters = global.cluster(params.threshold, params.grow);
    std::cout << "Found " << params.num_clusters << " clusters. Starting per-cluster analysis..." << std::endl;

    // Save the global classification volume once
    global.cluster_id.write(params.out_dir + "/" + params.filename + "_classes.bin");

    json report;
    report["filename"] = params.filename;
    report["dimension"] = Dimension;
    report["grid_size"] = params.grid;
    report["total_events"] = events.size();
    report["num_clusters"] = params.num_clusters;
    report["clusters"] = json::array();

    std::ofstream xml_out(params.filename + ".report.xml");
    xml_out << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
    xml_out << "<DarwinReport file=\"" << params.filename << "\" dimensions=\"" << Dimension << "\">\n";
    xml_out << "  <Summary totalEvents=\"" << events.size() << "\" clustersFound=\"" << params.num_clusters << "\" />\n";
    xml_out << "  <Clusters>\n";

    // Prepare HTML report
    std::ofstream html_out(params.out_dir + "/index.html");
    html_out << "<html><head><title>Darwin Report: " << params.filename << "</title>";
    html_out << "<style>body{font-family:sans-serif;margin:40px;} table{border-collapse:collapse;width:100%;} th,td{padding:10px;border:1px solid #ccc;text-align:left;} img{max-width:800px;}</style></head>";
    html_out << "<body><h1>Darwin Report: " << params.filename << "</h1>";
    html_out << "<table><tr><th>ID</th><th>Events</th><th>Diff Error</th><th>Factor Prob</th><th>Mean Vector</th><th>Visualization</th></tr>";

    std::string img_dir = params.out_dir + "/images";
    fs::create_directories(img_dir);

    std::vector<pid_t> active_pids;
    unsigned max_concurrent = std::max(1u, std::thread::hardware_concurrency());

    auto dispatch_plot = [&](const std::string &path, const std::vector<std::vector<double>> &class_data, const std::vector<std::vector<double>> &quant_data) {
        std::cout.flush();
        html_out.flush();
        xml_out.flush();

        while (active_pids.size() >= max_concurrent) {
            int status;
            pid_t done = waitpid(-1, &status, 0);
            if (done > 0) {
                active_pids.erase(std::remove(active_pids.begin(), active_pids.end(), done), active_pids.end());
            }
        }

        pid_t pid = fork();
        if (pid == 0) {
            make_plot(path, params, class_data, quant_data);
            std::exit(0);
        } else if (pid > 0) {
            active_pids.push_back(pid);
        }
    };

    std::vector<std::vector<typename Weighty<Dimension>::Event>> cluster_events(params.max_clusters + 1);
    for (const auto &e : events)
    {
        unsigned short c = 0;
        Coordinate<Dimension> coord(global);
        if (global.locate(e, coord))
        {
            size_t idx = (size_t)coord;
            c = static_cast<const Function<Dimension, unsigned short> &>(global.cluster_id)[idx];
        }
        if (c < params.max_clusters)
            cluster_events[c].push_back(e);
        else
            cluster_events[0].push_back(e);
    }
    for_each_plane<Dimension>([&params, &cluster_events, &marginal, &img_dir, &dispatch_plot](unsigned i, unsigned j)
    {
        marginal.reset();
        Coordinate<2> marginal_coord(marginal);
        Weighty<2>::Event marginal_event;
        std::vector<unsigned short> marginal_klass(params.points * params.points, 0);

        for (unsigned short c = 0; c <= params.num_clusters; ++c)
        {
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

            std::string path = img_dir + "/sample_" + std::to_string(i + 1) + "X" + std::to_string(j + 1) + ".png";
            dispatch_plot(path, class_data, quant_data);
    } });
    
    for (unsigned c = 1; c <= std::min(params.num_clusters, params.max_clusters); ++c)
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

        for_each_plane<Dimension>([&params, c, &cluster_events, &marginal, &img_dir, &dispatch_plot](unsigned i, unsigned j)
        {
            Weighty<2>::Event marginal_event;
            Coordinate<2> marginal_coord(marginal);
            std::vector<unsigned short> marginal_klass(params.points * params.points, 0);

            marginal.reset();
            for (auto &e : cluster_events[c])
            {
                marginal_event[0] = e[i];
                marginal_event[1] = e[j];
                marginal.event(marginal_event);
                if (marginal.locate(marginal_event, marginal_coord))
                {
                    size_t idx = (size_t)marginal_coord;
                    marginal_klass[idx] = c;
                }
            }
            marginal.prepare(params.smooth);

            // make cluster marginal plot
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

                std::string path = img_dir + "/cluster_" + std::to_string(c) + "_" + std::to_string(i + 1) + "X" + std::to_string(j + 1) + ".png";
                dispatch_plot(path, class_data, quant_data);
        } });

        std::string img_rel_name = "images/cluster_" + std::to_string(c) + "_1X2.png";

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
    }

    if (!active_pids.empty()) {
        std::cout << "Waiting for " << active_pids.size() << " background plots to finish rendering..." << std::endl;
        for (pid_t pid : active_pids) {
            int status;
            waitpid(pid, &status, 0);
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

    // Force Gnuplot to use a non-interactive terminal. 
    setenv("GNUTERM", "png", 1);

    cxxopts::Options options("Darwin", "Hierarchical Laplacian and Riemannian analysis");

    options.add_options()
    ("f,file", "Input filename", cxxopts::value<std::string>()->default_value("test_data"))
    ("d,dimension", "Dimension", cxxopts::value<unsigned>()->default_value("2"))
    ("g,grid", "Grid resolution", cxxopts::value<unsigned>())
    ("s,smooth", "Smoothing factor", cxxopts::value<float>()->default_value("0.01"))
    ("t,threshold", "Threshold", cxxopts::value<float>()->default_value("0.001"))
    ("visual", "Save visualization files", cxxopts::value<bool>()->default_value("false")->implicit_value("true"))
    ("max-clusters", "Maximum number of clusters to report", cxxopts::value<unsigned>()->default_value("10"))
    ("min-events", "Minimum number of events for a cluster to be reported", cxxopts::value<size_t>()->default_value("1"))
    ("v,verbose", "Verbose", cxxopts::value<bool>()->default_value("false")->implicit_value("true"))
    ("antialias", "Antialiasing", cxxopts::value<bool>()->default_value("true")->implicit_value("true"))
    ("verify", "Verify consistency", cxxopts::value<bool>()->default_value("true")->implicit_value("true"))
    ("grow", "Grow clusters", cxxopts::value<bool>()->default_value("false")->implicit_value("true"))
    ("a,ascii", "Use ASCII data", cxxopts::value<bool>()->default_value("false")->implicit_value("true"))
    ("h,help", "Print usage");

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
