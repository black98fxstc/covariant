#include <iostream>
#include <string>
#include <fstream>
#include <cxxopts.hpp>
#include <functional>

#include "Covariant.hpp"
#include "Samples.hpp"

struct Params
{
    std::string filename;
    std::vector<std::string> variables;
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
int do_it(const Params &params, unsigned grid_size, const Projection& proj)
{
    std::cout << "Laplace running with" 
              << " filename=" << params.filename << " smooth=" << params.smooth << " threshold=" << params.threshold
              << " dimension=" << params.dimension << " grid=" << params.grid << std::endl;
    std::cout << "   grow=" << (params.grow ? "on" : "off") << " verify=" << (params.verify ? "on" : "off") << " antialias=" << (params.antialias ? "on" : "off") 
              << " verbose=" << (params.verbose ? "on" : "off") << " visual=" << (params.visual ? "on" : "off") << std::endl;

    Events<Dimension> events;
    size_t num_events = proj[0].get_data<float>().size();
    events.resize(num_events);
    for (unsigned d = 0; d < Dimension; ++d) {
        const auto& col = proj[d].get_data<float>();
        for (size_t i = 0; i < num_events; ++i) {
            events[i][d] = col[i];
        }
    }
    std::cout << "Loaded " << events.size() << " events." << std::endl;

    std::cout << "Processing " << events.size() << " events..." << std::endl;
    Laplace<Dimension> laplace(grid_size);
    laplace.visualize = params.visual;
    laplace.antialias = params.antialias;
    laplace.verify = params.verify;
    laplace.verbose = params.verbose;
    for (const auto &e : events)
        laplace.event(e);
    std::cout << "Found " << events.size() << " valid events..." << std::endl;

    std::cout << "Analyzing the sample..." << std::endl;
    laplace.analyze(params.smooth, params.threshold);
    if (laplace.differentialError() > .0001)
        std::cout << "Differential equation solution is unusually bad " << laplace.differentialError() << std::endl;
    else
        std::cout << "Consistency checks passed..." << std::endl;

    std::cout << "Performing Laplacian clustering..." << std::endl;
    unsigned found = laplace.cluster(params.threshold, params.grow);
    std::vector<unsigned short> classes;
    classes.reserve(events.size());
    for (const auto &e : events)
        classes.push_back(laplace.classify(e));
    std::cout << "Found " << found << " clusters." << std::endl;

    if (params.ascii)
    {
        std::ofstream out(params.filename + ".cluster");
        if (!out)
        {
            std::cerr << "Error: Could not open cluster file for writing: " << params.filename + ".cluster" << std::endl;
            return 1;
        }
        for (auto c : classes)
            out << c << "\n";
        out.close();
    }
    else
    {
        std::ofstream out(params.filename + ".cluster", std::ios::binary | std::ios::trunc);
        if (!out)
        {
            std::cerr << "Error: Could not open cluster file for writing: " << params.filename + ".cluster" << std::endl;
            return 1;
        }
        out.write(reinterpret_cast<const char *>(classes.data()), classes.size() * sizeof(unsigned short));
        out.close();
    }
    std::cout << "Saved to file " << params.filename + ".cluster" << std::endl;

    return 0;
}

int main(int argc, char *argv[])
{
    Params params;

    // Set up the command-line options parser.
    cxxopts::Options options("CovariantCLI", "A command-line interface for the Laplacian clustering algorithm");

    // Add the command-line options.
    options.add_options()("f,file", "File name for data", cxxopts::value<std::string>()->default_value("test_data"))
    ("variables", "List of variables", cxxopts::value<std::vector<std::string>>())
    ("populations", "List of populations", cxxopts::value<std::vector<std::string>>())
    ("d,dimension", "Dimension of the events", cxxopts::value<unsigned>()->default_value("2"))
    ("g,grid", "Grid resolution", cxxopts::value<unsigned>())
    ("s,smooth", "Smoothing factor", cxxopts::value<float>()->default_value("0.01"))
    ("t,threshold", "Consistency threshold", cxxopts::value<float>()->default_value("0.001"))
    ("visual", "Save files for visualization", cxxopts::value<bool>()->default_value("false")->implicit_value("true"))
    ("v,verbose", "Verbose output", cxxopts::value<bool>()->default_value("false")->implicit_value("true"))
    ("antialias", "Enable antialiasing (developer feature)", cxxopts::value<bool>()->default_value("true"))
    ("verify", "Enable consistency verification", cxxopts::value<bool>()->default_value("true")->implicit_value("true"))
    ("grow", "Unclaimed events are assigned to the nearest cluster", cxxopts::value<bool>()->default_value("false")->implicit_value("true"))
    ("a,ascii", "Use ASCII format for data files", cxxopts::value<bool>()->default_value("false")->implicit_value("true"))
    ("h,help", "Print usage");

    options.parse_positional({"file", "variables", "populations"});

    auto result = options.parse(argc, argv);

    if (result.count("help"))
    {
        std::cout << options.help() << std::endl;
        return 0;
    }

    params.filename = result["file"].as<std::string>();
    if (result.count("variables")) {
        params.variables = result["variables"].as<std::vector<std::string>>();
        params.dimension = params.variables.size();
    } else {
        params.dimension = result["dimension"].as<unsigned>();
    }
    params.smooth = result["smooth"].as<float>();
    params.threshold = result["threshold"].as<float>();
    params.visual = result["visual"].as<bool>();
    params.antialias = result["antialias"].as<bool>();
    params.verify = result["verify"].as<bool>();
    params.grow = result["grow"].as<bool>();
    params.verbose = result["verbose"].as<bool>();
    params.ascii = result["ascii"].as<bool>();

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
        case 4:
            params.grid = 32;
            break;
        default:
            params.grid = 32;
            break;
        }

    DataSet dataset;
    if (!dataset.read(params.filename)) {
        std::cerr << "Error: Could not open data file: " << params.filename << std::endl;
        return 1;
    }

    std::vector<size_t> var_indices;
    if (!params.variables.empty()) {
        for (const auto& var_name : params.variables) {
            bool found = false;
            for (size_t i = 0; i < dataset.size(); ++i) {
                if (dataset[i].name == var_name) {
                    var_indices.push_back(i);
                    found = true;
                    break;
                }
            }
            if (!found) {
                std::cerr << "Error: Could not find variable '" << var_name << "' in dataset." << std::endl;
                return 1;
            }
        }
    } else {
        for (size_t i = 0; i < params.dimension && i < dataset.size(); ++i) {
            var_indices.push_back(i);
        }
        if (var_indices.size() != params.dimension) {
            std::cerr << "Error: Dataset has fewer variables than requested dimension." << std::endl;
            return 1;
        }
    }

    Projection proj(var_indices, dataset);

    std::vector<size_t> pop_indices;
    if (result.count("populations")) {
        for (const auto& pop_name : result["populations"].as<std::vector<std::string>>()) {
            bool found = false;
            for (size_t i = 0; i < dataset.size(); ++i) {
                if (dataset[i].name == pop_name) {
                    pop_indices.push_back(i);
                    found = true;
                    break;
                }
            }
            if (!found) {
                std::cerr << "Error: Could not find population '" << pop_name << "' in dataset." << std::endl;
                return 1;
            }
        }
    }

    switch (params.dimension)
    {
    case 2:
        return do_it<2>(params, params.grid, proj);
    case 3:
        return do_it<3>(params, params.grid, proj);
    case 4:
        return do_it<4>(params, params.grid, proj);
    default:
        std::cerr << "Error: Unsupported dimension: " << params.dimension << std::endl;
        return 1;
    }
};