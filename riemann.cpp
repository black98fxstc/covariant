#include <iostream>
#include <string>
#include <fstream>
#include <functional>
#include <cxxopts.hpp>

#include "Covariant.hpp"
#include "Samples.hpp"

struct Params
{
    std::string filename;
    std::vector<std::string> variables;
    unsigned dimension;
    unsigned cluster;
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
int do_it(const Params &params, const cxxopts::ParseResult &result, unsigned grid_size, const Projection& proj)
{
    std::cout << "Riemann running with" 
              << " filename=" << params.filename << " smooth=" << params.smooth << " threshold=" << params.threshold
              << " dimension=" << params.dimension << " grid=" << params.grid << std::endl;
    std::cout << "   grow=" << (params.grow ? "on" : "off") << " verify=" << (params.verify ? "on" : "off") << " antialias=" << (params.antialias ? "on" : "off") 
              << " verbose=" << (params.verbose ? "on" : "off") << " visual=" << (params.visual ? "on" : "off") << std::endl;

    typename Weighty<Dimension>::Events events;
    size_t num_events = proj[0].get_data<float>().size();
    events.resize(num_events);
    for (unsigned d = 0; d < Dimension; ++d) {
        const auto& col = proj[d].get_data<float>();
        for (size_t i = 0; i < num_events; ++i) {
            events[i][d] = col[i];
        }
    }
    std::cout << "Loaded " << events.size() << " events." << std::endl;

    Riemann<Dimension> riemann(grid_size); // Use passed grid_size
    riemann.visualize = params.visual;
    riemann.antialias = params.antialias;
    riemann.verify = params.verify;
    riemann.verbose = params.verbose;
    size_t valid_events = 0;
    if (result.count("cluster"))
    {
        std::vector<unsigned short> classes(events.size());
        if (params.ascii)
        {
            std::ifstream in(params.filename + ".cluster");
            for (unsigned i = 0; i < events.size(); i++)
                in >> classes[i];
            in.close();
        }
        else
        {
            std::ifstream in(params.filename + ".cluster", std::ios::binary);
            in.read(reinterpret_cast<char *>(classes.data()), events.size() * sizeof(unsigned short));
            in.close();
        }

        std::cout << "Processing cluster " << params.cluster << " of " << events.size() << " events..." << std::endl;
        for (unsigned i = 0; i < events.size(); i++)
        {
            if (classes[i] != params.cluster)
                continue;
            if (riemann.event(events[i]))
                valid_events++;
        }
    }
    else
    {
        std::cout << "Processing " << events.size() << " events..." << std::endl;
        for (const auto &e : events)
            if (riemann.event(e))
                valid_events++;
    }
    std::cout << "Found " << events.size() << " valid events..." << std::endl;

    std::cout << "Analyzing the sample..." << std::endl;

    riemann.analyze(params.smooth, params.threshold);
    if (riemann.factorProbability() > .0001)
    {
        std::cout << "Probability factoring is unusually bad " << riemann.factorProbability() << std::endl; // No change here, it's a method call
    }
    else if (riemann.differentialError() > .0001)
    {
        std::cout << "Differential equation solution is unusually bad " << riemann.differentialError() << std::endl;
    }
    else
    {
        std::cout << "Consistency checks passed..." << std::endl;
    }

    if (!result.count("cluster"))
    {
        std::cout << "Performing Laplacian clustering..." << std::endl;
        unsigned found = riemann.cluster(params.threshold, params.grow);
        std::vector<unsigned short> classes;
        classes.reserve(events.size());
        for (const auto &e : events)
            classes.push_back(riemann.classify(e));
        std::cout << "Found " << found << " clusters..." << std::endl;

        if (params.ascii)
        {
            std::ofstream out(params.filename + ".cluster");
            for (auto c : classes)
                out << c << "\n";
            out.close();
        }
        else
        {
            std::ofstream out(params.filename + ".cluster", std::ios::binary | std::ios::trunc);
            out.write(reinterpret_cast<const char *>(classes.data()), classes.size() * sizeof(unsigned short));
            out.close();
        }
        std::cout << "Saved to file " << params.filename + ".cluster" << std::endl;
    }
    return 0; // Return 0 on success
}

int main(int argc, char *argv[])
{
    Params params;

    // Set up the command-line options parser.
    cxxopts::Options options("CovariantCLI", "Generate test data for the Covariant class");

    // Add the command-line options.
    options.add_options()
    ("f,file", "Output filename for generated data", cxxopts::value<std::string>()->default_value("test_data"))
    ("variables", "List of variables", cxxopts::value<std::vector<std::string>>())
    ("populations", "List of populations", cxxopts::value<std::vector<std::string>>())
    ("c,cluster", "Process data from the specified cluster", cxxopts::value<unsigned>()->default_value("0"))
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
    params.filename = result["file"].as<std::string>();
    params.cluster = result["cluster"].as<unsigned>();

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
        return do_it<2>(params, result, params.grid, proj);
    case 3:
        return do_it<3>(params, result, params.grid, proj);
    case 4:
        return do_it<4>(params, result, params.grid, proj);

    default:
        std::cerr << "Error: Unsupported dimension: " << params.dimension << std::endl;
        return 1;
    }
    return 0; // Should not be reached if all cases return
}