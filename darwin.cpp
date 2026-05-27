#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <cxxopts.hpp>
#include <nlohmann/json.hpp>
#include <cstdlib>

#include "Covariant.hpp"

using json = nlohmann::json;

struct Params
{
    std::string filename;
    unsigned dimension;
    unsigned grid;
    float smooth;
    float threshold;
    bool visual;
    bool verbose;
    bool antialias;
    bool verify;
    bool grow;
    bool ascii;
};

template <unsigned Dimension>
int do_it(const Params &params, unsigned grid_size)
{
    TestData<Dimension> events;
    std::string ext = params.ascii ? ".txt" : ".dat";
    std::cout << "Loading events from " << params.filename << ext << "..." << std::endl;
    if (!events.read(params.filename + ext, params.ascii))
    {
        std::cerr << "Error: Could not open event file." << std::endl;
        return 1;
    }

    // 1. Initial Riemann object to perform global clustering via Laplace
    Riemann<Dimension> global(grid_size, true);
    global.visualize = params.visual;
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

    json report;
    report["filename"] = params.filename;
    report["dimension"] = Dimension;
    report["grid_size"] = grid_size;
    report["total_events"] = events.size();
    report["num_clusters"] = num_clusters;
    report["clusters"] = json::array();

    std::ofstream xml_out(params.filename + ".report.xml");
    xml_out << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
    xml_out << "<DarwinReport file=\"" << params.filename << "\" dimensions=\"" << Dimension << "\">\n";
    xml_out << "  <Summary totalEvents=\"" << events.size() << "\" clustersFound=\"" << num_clusters << "\" />\n";
    xml_out << "  <Clusters>\n";

    for (unsigned c = 1; c <= num_clusters; ++c)
    {
        std::cout << "  Analyzing cluster " << c << "..." << std::endl;
        size_t cluster_event_count = 0;
        global.reset();
        for (const auto &e : events)
        {
            if (global.classify(e) == c)
            {
                global.event(e);
                cluster_event_count++;
            }
        }

        global.Riemann<Dimension>::analyze(params.smooth, params.threshold);

        // Populate JSON data
        json c_info;
        c_info["id"] = c;
        c_info["event_count"] = cluster_event_count;
        c_info["differential_error"] = global.differentialError();
        c_info["factor_probability"] = global.factorProbability();
        report["clusters"].push_back(c_info);

        std::string img_name = params.filename + ".cluster" + std::to_string(c) + ".png";

        // Populate XML data
        xml_out << "    <Cluster id=\"" << c << "\" events=\"" << cluster_event_count << "\">\n";
        xml_out << "      <Metrics diffError=\"" << global.differentialError() << "\" ";
        xml_out << "factorProb=\"" << global.factorProbability() << "\" />\n";
        xml_out << "      <Visualization src=\"" << img_name << "\" />\n";
        xml_out << "    </Cluster>\n";

        if (params.visual)
        {
            std::cout << "      Generating visualization: " << img_name << "..." << std::endl;
            std::string cmd = "matlab -batch \"run('covariant2.m'); exportgraphics(gcf, '" + img_name + "', 'Resolution', 150);\"";
            std::system(cmd.c_str());
        }
    }

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

    options.add_options()("f,file", "Input filename", cxxopts::value<std::string>()->default_value("test_data"))("d,dimension", "Dimension", cxxopts::value<unsigned>()->default_value("2"))("g,grid", "Grid resolution", cxxopts::value<unsigned>())("s,smooth", "Smoothing factor", cxxopts::value<float>()->default_value("0.01"))("t,threshold", "Threshold", cxxopts::value<float>()->default_value("0.001"))("visual", "Save visualization files", cxxopts::value<bool>()->default_value("false"))("v,verbose", "Verbose", cxxopts::value<bool>()->default_value("false"))("antialias", "Antialiasing", cxxopts::value<bool>()->default_value("true"))("verify", "Verify consistency", cxxopts::value<bool>()->default_value("true"))("grow", "Grow clusters", cxxopts::value<bool>()->default_value("false"))("a,ascii", "Use ASCII data", cxxopts::value<bool>()->default_value("false"))("h,help", "Print usage");

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
    params.antialias = result["antialias"].as<bool>();
    params.verify = result["verify"].as<bool>();
    params.grow = result["grow"].as<bool>();
    params.verbose = result["verbose"].as<bool>();
    params.ascii = result["ascii"].as<bool>();

    if (result.count("grid"))
    {
        params.grid = result["grid"].as<unsigned>();
    }
    else
    {
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
    }

    std::cout << "Darwin running: " << params.filename << " (Dim: " << params.dimension << ")" << std::endl;

    switch (params.dimension)
    {
    case 2:
        return do_it<2>(params, params.grid);
    case 3:
        return do_it<3>(params, params.grid);
    case 4:
        return do_it<4>(params, params.grid);
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
