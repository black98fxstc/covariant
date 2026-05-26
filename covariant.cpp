#include <iostream>
#include <string>
#include <fstream>
#include <functional>
#include <cxxopts.hpp>

#include "Covariant.hpp"

struct Params
{
    std::string filename;
    unsigned dimension;
    unsigned cluster;
    unsigned grid;
    float smooth;
    float threshold;
    bool visual;
    bool verbose;
    bool antialias;
    bool verify;
    bool ascii;
};

template <unsigned Dimension>
int do_it(const Params &params, const cxxopts::ParseResult &result, unsigned grid_size)
{
    TestData<Dimension> events;
    std::string ext = params.ascii ? ".txt" : ".dat";
    std::cout << "Loading events from " << params.filename << ext << "..." << std::endl;
    if (!events.read(params.filename + ext, params.ascii))
    {
        std::cerr << "Error: Could not open event file for loading: " << params.filename + ext << std::endl;
        return 1; // Return 1 on error
    }
    std::cout << "Loaded " << events.size() << " events." << std::endl;

    Covariant<Dimension> covariant(grid_size, true); // Use passed grid_size
    covariant.visualize = params.visual;
    covariant.antialias = params.antialias;
    covariant.verify = params.verify;
    covariant.verbose = params.verbose;
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
            if (covariant.event(events[i]))
                valid_events++;
        }
    }
    else
    {
        std::cout << "Processing " << events.size() << " events..." << std::endl;
        for (const auto &e : events)
            if (covariant.event(e))
                valid_events++;
    }
    std::cout << "Found " << events.size() << " valid events..." << std::endl;

    std::cout << "Analyzing the sample..." << std::endl;

    covariant.analyze(params.smooth, params.threshold);
    if (covariant.factorProbability() > .0001)
    {
        std::cout << "Probability factoring is unusually bad " << covariant.factorProbability() << std::endl; // No change here, it's a method call
    }
    else if (covariant.differentialEquation() > .0001)
    {
        std::cout << "Differential equation solution is unusually bad " << covariant.differentialEquation() << std::endl;
    }
    else
    {
        std::cout << "Consistency checkes passed..." << std::endl;
    }

    if (!result.count("cluster"))
    {
        std::cout << "Performing Laplacian clustering..." << std::endl;
        unsigned found = covariant.cluster(params.threshold);
        std::vector<unsigned short> classes;
        classes.reserve(events.size());
        for (const auto &e : events)
            classes.push_back(covariant.classify(e));
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
    options.add_options()("f,file", "Output filename for generated data", cxxopts::value<std::string>()->default_value("test_data"))("c,cluster", "Process data from the specified cluster", cxxopts::value<unsigned>()->default_value("0"))("d,dimension", "Dimension of the events", cxxopts::value<unsigned>()->default_value("2"))("g,grid", "Grid resolution", cxxopts::value<unsigned>())("s,smooth", "Smoothing factor", cxxopts::value<float>()->default_value("0.01"))("t,threshold", "Consistency threshold", cxxopts::value<float>()->default_value("0.001"))("visual", "Enable visualization", cxxopts::value<bool>()->default_value("false"))("v,verbose", "Verbose output", cxxopts::value<bool>()->default_value("false"))("antialias", "Enable antialiasing (developer feature)", cxxopts::value<bool>()->default_value("true"))("verify", "Enable consistency verification", cxxopts::value<bool>()->default_value("true"))("a,ascii", "Use ASCII format for data files", cxxopts::value<bool>()->default_value("false"))("h,help", "Print usage");

    options.parse_positional({"file"});

    auto result = options.parse(argc, argv);

    if (result.count("help"))
    {
        std::cout << options.help() << std::endl;
        return 0;
    }

    params.filename = result["file"].as<std::string>();
    params.dimension = result["dimension"].as<unsigned>();
    params.cluster = result["cluster"].as<unsigned>();

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
        case 4:
            params.grid = 32;
            break;
        default:
            params.grid = 32;
            break;
        }
    }

    params.smooth = result["smooth"].as<float>();
    params.threshold = result["threshold"].as<float>();
    params.visual = result["visual"].as<bool>();
    params.antialias = result["antialias"].as<bool>();
    params.verify = result["verify"].as<bool>();
    params.verbose = result["verbose"].as<bool>();
    params.ascii = result["ascii"].as<bool>();

    std::cout << "Program running with dimension=" << params.dimension << " grid=" << params.grid << " filename=" << params.filename
              << " smooth=" << params.smooth << " threshold=" << params.threshold
              << " antialias=" << (params.antialias ? "on" : "off")
              << " verify=" << (params.verify ? "on" : "off") << std::endl;
    std::cout << "Analyzing events in " << params.dimension << " dimensions..." << std::endl;

    switch (params.dimension)
    {
    case 2:
        return do_it<2>(params, result, params.grid);
    case 3:
        return do_it<3>(params, result, params.grid);
    case 4:
        return do_it<4>(params, result, params.grid);

    default:
        std::cerr << "Error: Unsupported dimension: " << params.dimension << std::endl;
        return 1;
    }
    return 0; // Should not be reached if all cases return
}